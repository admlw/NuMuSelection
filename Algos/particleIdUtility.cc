#include "particleIdUtility.h"

namespace pidutil{

  std::pair<double, double> particleIdUtility::getAveragedQdX(recob::Track const& track, std::vector< art::Ptr< recob::Hit > > hitCollection, std::vector< art::Ptr< raw::RawDigit > > rawDCollection, bool isDqdxFromRawD){

    std::pair<double, double> averagedQdX;
    averagedQdX.first = -1;

    std::pair<int, int> channelPeakTime;
    std::vector< std::pair< int, int > > channelPeakTimeVector;

    // get total charge deposited
    double hitCharge = 0.0;
    double rawIntegral = 0.0;
    std::vector<double> hitCharges;
    for (size_t i = 0; i < hitCollection.size(); i++){

      art::Ptr<recob::Hit> hit = hitCollection.at(i);

      // select out "good" hits on the Y plane
      if (hit->Channel() < 4800 || hit->GoodnessOfFit() > 1.1) continue;

      channelPeakTime.first = hit->Channel();
      channelPeakTime.second = hit->PeakTime();
      channelPeakTimeVector.push_back(channelPeakTime);

      // get RawDigitIntegral
      for (size_t j = 0; j < rawDCollection.size(); j++){

        if ((int)rawDCollection.at(j)->Channel() == (int)channelPeakTime.first){

          // get correction factor
          std::vector<short> const adcVec = rawDCollection.at(j)->ADCs();
          double correctionValue = 0;
          int counter = 0;
          for ( int k = 0; k < (int)adcVec.size(); k++){

            if ( k < channelPeakTime.second - 30 || k > channelPeakTime.second + 30 ){

              counter++;
              correctionValue +=adcVec.at(k);

            }

          }

          correctionValue = correctionValue/(double)counter;

          for ( int k = channelPeakTime.second - 30; 
              k < channelPeakTime.second + 30;
              k++){

            rawIntegral += (rawDCollection.at(j)->ADC(k)) - correctionValue;

          }
        }
      }


      hitCharge += hit->Integral();
      hitCharges.push_back(hit->Integral());

    }

    std::sort(channelPeakTimeVector.begin(), channelPeakTimeVector.end(), pairCompare);

    std::vector< std::vector< std::pair<double, double> > > deadRegions =
      getDeadRegions(channelPeakTimeVector);

    double untrackedLength = getUntrackedLength(track, deadRegions);

    // get full track length
    double trackLength = track.Length() - untrackedLength;


    double averagedQdXmean; 
    if (isDqdxFromRawD == false){
      averagedQdXmean = hitCharge/trackLength;
    }
    else {
      averagedQdXmean = rawIntegral/trackLength;
    }

    averagedQdX.first = averagedQdXmean;
    averagedQdX.second = trackLength;

    return averagedQdX;

  }

  std::vector< std::vector< std::pair<double, double> > > particleIdUtility::getDeadRegions(std::vector< std::pair<int, int> > hitCollection){

    const detinfo::DetectorProperties* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    std::pair<double, double> deadStart;
    std::pair<double, double> deadEnd;
    std::vector< std::pair<double, double> > deadStartVec;
    std::vector< std::pair<double, double> > deadEndVec;

    for (size_t i = 1; i < hitCollection.size(); i++){

      std::pair<int,int> prevHit = hitCollection.at(i-1);
      std::pair<int,int> thisHit = hitCollection.at(i);

      if (std::abs(prevHit.first - thisHit.first) > 1){

        deadStart.first = (prevHit.first-4800)*0.3;
        deadStart.second = detProp->ConvertTicksToX(prevHit.second, geo::PlaneID(0,0,2));
        deadStartVec.push_back(deadStart);

        deadEnd.first = (thisHit.first-4800)*0.3;
        deadEnd.second = detProp->ConvertTicksToX(prevHit.second,geo::PlaneID(0,0,2));
        deadEndVec.push_back(deadEnd);

      }

    }

    std::vector< std::vector< std::pair< double, double > > > returnVec;
    returnVec.push_back(deadStartVec);
    returnVec.push_back(deadEndVec);

    return returnVec;

  }

  double particleIdUtility::getUntrackedLength(recob::Track const& track, std::vector< std::vector< std::pair<double,double> > > deadRegions){

    std::vector< std::pair<double,double> > deadStart = deadRegions.at(0);
    std::vector< std::pair<double,double> > deadEnd = deadRegions.at(1);

    double untrackedLength = 0.0;

    // loop over number of dead regions
    for ( size_t i = 0; i < deadStart.size(); i++){


      int deadStartPoint = -1;
      int deadEndPoint = -1;

      // get particle trajectory points
      for (size_t j = 0; j < track.NumberTrajectoryPoints(); j++){

        Point_t trajPoint = track.TrajectoryPoint(j).position;

        if( trajPoint.Z() >= deadStart.at(i).first && deadStartPoint == -1){

          deadStartPoint = j;

        }

        if ( trajPoint.Z() <= deadEnd.at(i).first && deadStartPoint != -1){

          deadEndPoint = j;

        }

      }

      std::cout << "Found dead region between trajectory points " 
        << deadStartPoint << " and " << deadEndPoint 
        << " which corresponds to z positions " << deadStart.at(i).first 
        << " and " << deadEnd.at(i).first << std::endl;

      // calculate length from trajectory in dead region
      if (deadStartPoint >= deadEndPoint) return 0;

      size_t iCurr = track.Trajectory().NextValidPoint((size_t)deadStartPoint);
      size_t iNext = iCurr;
      size_t iLast = (size_t)deadEndPoint;
      Point_t const* curr = &(track.Trajectory().LocationAtPoint(iCurr));
      Coord_t len = 0.0;
      while ((iNext = track.Trajectory().NextValidPoint(++iNext)) <= iLast){

        Point_t const* next = &track.Trajectory().LocationAtPoint(iNext);
        len += (*next - *curr).R();
        curr = next;

      }
      untrackedLength += len; 

    }

    std::cout << "Returning untracked Length " << untrackedLength << std::endl;

    return untrackedLength;

  }

  bool particleIdUtility::isMuonLike(double averagedQdX, double len){

    if (averagedQdX > 300){
      TF1* f1 = new TF1("f1", "500000/(100*(x-300))-10", 300, 2000);

      double evalPoint = f1->Eval(averagedQdX);

      std::cout << "\n average dqdx " << averagedQdX
        << "\n length " << len 
        << "\n evalPoint " << evalPoint
        << "\n " << std::endl;

      if (len > evalPoint) return false;
      else return true;
    }
    else return true;
  }
}
