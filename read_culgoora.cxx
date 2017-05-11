////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  analysis.cxx
////  Analyzing Solar Flare Data
////
////  March 2017
////  Want to make dynamic spectrogram
////  Trying to use the ROOT time axis method like Amy wants
////  Now trying to do background subtraction to make the spectrogram look cleaner
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <memory>
#include <boost/scoped_ptr.hpp>
#include <vector>
//#include <array>
#include <cmath>
#include <algorithm>
//#include <complex.h>
#include <fftw3.h>
#include <ctime>
#include "time.h" // for time convert
#include <iterator>
#include <complex>

//for culgoora analysis
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

//AraRoot Includes
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "RawAraStationEvent.h"
#include "AraEventCorrelator.h"
#include "AraAntennaInfo.h"
#include "AraGeomTool.h"
#include "FFTtools.h"
#include "AraStationInfo.h"
#include "AraAntennaInfo.h"

//ROOT Includes
#include "TTreeIndex.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TText.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"
#include "TString.h"
#include "TAttText.h"
#include "TRandom.h"
#include "TTimeStamp.h"

//removing lots of auxilliary functions, I'm tired of having them clog up my code...

TGraph *makeSpectrum_mVPerRootHz(TGraph *grWave);

using namespace std;

//AraGeomTool *geomTool;
//RawIcrrStationEvent *rawIcrrEvPtr=0;
//RawAtriStationEvent *rawAtriEvPtr=0;
//RawAraStationEvent *rawEvPtr=0;
//UsefulIcrrStationEvent *realIcrrEvPtr=0;
//UsefulAtriStationEvent *realAtriEvPtr=0;

/* Channel mappings for the testbed
Channel 0: H Pol
Channel 1: H Pol
Channel 2: V Pol
Channel 3: V Pol
Channel 4: V Pol
Channel 5: H Pol
Channel 6: V Pol
Channel 7: H Pol
Channel 8: V Pol
Channel 9: H Pol
Channel 10: V Pol
Channel 11: H Pol
Channel 12: H Pol
Channel 13: H Pol
Channel 14: Surface
Channel 15: Surface
*/
 
 

int main(int argc, char **argv)
{
	int fd;
	char c;
	
	fd = open("SPEC110215", O_RDONLY);
	while (read(fd, &c, 1) > 0) {
		printf("%d\n", c);
	}
	close(fd);
	
	
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
    
	std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <input file 1> <input file 2> ....<input file n>" << std::endl;
		return -1;
	}
	
	gStyle->SetOptStat(0000000000);

	time_t time_now = time(0); //get the time now
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;
		
	int num_found =0; //a variable to hold the number of events we've found in the first file
	int num_found_max=100; //a variable to hold the maximum number of waveforms you want to allow in the first input file loop
	
	//TString dataPath("/home/brianclark/SolarFlares/"); //the location of the files
	TString dataPath("/n/home00/clark.2668/workspace/TrunkSolarFlares/source_AraRoot_Trunk/");
	TChain chain("eventTree"); //this for the events for the exterior loop
	for(int file=1; file<argc; file++){
		TString fileKey(argv[file]); //a file key
		chain.Add(dataPath+fileKey); //add files to the chain
	}
	
	RawIcrrStationEvent *rawIcrrEvPtr=0; //make the raw event pointer
	int RunNo=0;
	chain.SetBranchAddress("event",&rawIcrrEvPtr); //set the branch address
	chain.SetBranchAddress("run",&RunNo);
	Int_t numEntries = chain.GetEntries(); //get the number of entries
	cout<<"Num of entries "<<numEntries<<endl;
	//numEntries = 132;
	
	//okay, need to set up the TH2D for the final spectrogram
	//x axis is time
	
	TDatime normal_start(2011, 02, 15, 1, 00,0);
	int start_time_bin = normal_start.Convert()-18000;
	TDatime normal_stop(2011, 02, 15, 2, 55,0);
	int stop_time_bin = normal_stop.Convert()-18000;
	
	TDatime background_start(2011, 02, 11, 1, 00,0);
	int start_time_bin_background = background_start.Convert()-18000;
	TDatime background_stop(2011, 02, 11, 2, 55,0);
	int stop_time_bin_background = background_stop.Convert()-18000;
		
	/*
	int start_time_bin = 1297731600; 
	int stop_time_bin = start_time_bin + 10500;
	int num_time_bins = 3375;
	//gStyle -> SetTimeOffset(start_time_bin+14400-3600);
	
	int start_time_bin_background = 1297386000;
	int stop_time_bin_background = start_time_bin + 10500;
	//gStyle -> SetTimeOffset(start_time_bin_background+14400-3600);
	*/
	int num_time_bins = 3375;
	
	//y axis is frequency in MHz, binned in ~3.9 MHz bins
	//we'll interpolate to 0.5 ns intervals, which will give me 1 GHz of frequency bins to work with
	
	int num_freq_bins = 256;
	int start_freq_bin = 0;
	int stop_freq_bin = 1000;
	
	TH2D* spectrogram[16]; //a spectrogram to hold the Feb 15 events
	TH2D* spectrogram_background[16]; //a spectrogram to hold the Feb 11 background
	TH2D* spectrogram_cleaned[16]; //the difference between the two
		
	for(int i=0; i<16; i++){
		
		//these are referenced to Feb 15
		spectrogram[i] = new TH2D("","",num_time_bins,start_time_bin,stop_time_bin,num_freq_bins,start_freq_bin,stop_freq_bin); //the spectrogram 
		spectrogram[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		spectrogram[i]->GetXaxis()->SetTimeFormat("%H:%M");
		//spectrogram[i]->GetXaxis()->SetTimeOffset(start_time_bin+14400-3600);
		
		spectrogram_cleaned[i] = new TH2D("","",num_time_bins,start_time_bin,stop_time_bin,num_freq_bins,start_freq_bin,stop_freq_bin); //the spectrogram 
		spectrogram_cleaned[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		spectrogram_cleaned[i]->GetXaxis()->SetTimeFormat("%H:%M");
		//spectrogram_cleaned[i]->GetXaxis()->SetTimeOffset(start_time_bin+14400-3600);
		
		//these are referenced to Feb 11
		spectrogram_background[i] = new TH2D("","",num_time_bins,start_time_bin_background,stop_time_bin_background,num_freq_bins,start_freq_bin,stop_freq_bin); //the spectrogram 
		spectrogram_background[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		spectrogram_background[i]->GetXaxis()->SetTimeFormat("%H:%M");
		//spectrogram[i]->GetXaxis()->SetTimeOffset(start_time_bin_background+14400-3600);

	}
	
	TH1D* number_of_entries  = new TH1D("","",num_time_bins, start_time_bin, stop_time_bin);
	TH1D* number_of_entries_background  = new TH1D("","",num_time_bins, start_time_bin_background, stop_time_bin_background);
	number_of_entries->GetXaxis()->SetTimeDisplay(1);
	number_of_entries_background->GetXaxis()->SetTimeDisplay(1);
	
	int num_found_background =0;
	int num_found_regular =0;
		
	for(int event=0;event<numEntries;event++) { //loop over all events
		
		//if(num_found>num_found_max) {continue;} //bounce out if you have found more than this many events
		chain.GetEvent(event);  //This line gets the RawIcrr or RawAtri Event
		
		int unixtime_1;
		unixtime_1 = rawIcrrEvPtr->head.unixTime; // get  event unixTime
		time_t test_time_1 = unixtime_1;
		tm *test_time_tm_1 = gmtime( &test_time_1 );
		int year = test_time_tm_1->tm_year+1900;
		int month = test_time_tm_1->tm_mon+1;
		int day = test_time_tm_1->tm_mday;
		int hour = test_time_tm_1->tm_hour;
		int min = test_time_tm_1->tm_min;
		int sec = test_time_tm_1->tm_sec;
				
		if( (hour >= 1) && (hour<3)){ //look for events in the 1-3 o'clock hour only
			
			UsefulIcrrStationEvent *realIcrrEvPtr = new UsefulIcrrStationEvent(rawIcrrEvPtr, AraCalType::kLatestCalib); //make a real event pointer for the Icrr station
			
			//then sort by day
			if(day == 11){//it's a background event from Feb 11
				if( realIcrrEvPtr->isCalPulserEvent() != true && rawIcrrEvPtr -> trig.trigType != 68) {//check to see if it is NOT a cal pulser, and if it IS a software
					number_of_entries_background->Fill(unixtime_1);
					for(int chan=0; chan<16; chan++){ //loop over channels
						TGraph *waveform = realIcrrEvPtr->getGraphFromRFChan(chan); //get the waveform
						TGraph *waveform_interpolated =FFTtools::getInterpolatedGraph(waveform,0.5); //interpolate the waveform to 0.5 ns
						TGraph *spectrum = makeSpectrum_mVPerRootHz(waveform_interpolated); //make the spectrum of the waveform
						int length=spectrum->GetN(); //get the length of the spectrum
						for (int i=0; i<length; i++){ //loop over the spectrum
							double freq =0.;
							double mag =0.;
							spectrum->GetPoint(i,freq,mag); //for this point, get it's frequency, and magnitude
							spectrogram_background[chan]->Fill(unixtime_1,freq,((pow(mag/1000.,2)))); //Fill that bin; I think it should work to just drop in the magnitude as a weight, *I think*
						}
					}
					num_found_background++;
				}
			}
			else if(day == 15){//it's a flare event from Feb 15
			     	if( realIcrrEvPtr->isCalPulserEvent() != true && rawIcrrEvPtr -> trig.trigType != 68) {//check to see if it is NOT a cal pulser, and if it IS NOT a software
					number_of_entries->Fill(unixtime_1);
					for(int chan=0; chan<16; chan++){ //loop over channels
						TGraph *waveform = realIcrrEvPtr->getGraphFromRFChan(chan); //get the waveform
						TGraph *waveform_interpolated =FFTtools::getInterpolatedGraph(waveform,0.5); //interpolate the waveform to 0.5 ns
						TGraph *spectrum = makeSpectrum_mVPerRootHz(waveform_interpolated); //make the spectrum of the waveform
						int length=spectrum->GetN(); //get the length of the spectrum
						for (int i=0; i<length; i++){ //loop over the spectrum
							double freq =0.;
							double mag =0.;
							spectrum->GetPoint(i,freq,mag); //for this point, get it's frequency, and magnitude
							spectrogram[chan]->Fill(unixtime_1,freq,((pow(mag/1000.,2)))); //Fill that bin; I think it should work to just drop in the magnitude as a weight, *I think*
						}
					}
					num_found_regular++;
				}
			}
			delete realIcrrEvPtr; //delete this
		}
	}//close the outer loop over events
	
	cout<<"I found "<<num_found_regular<<" flare events" <<endl;
	cout<<"I found "<<num_found_background<<" background events" <<endl;
	
	//now, we have to normalize the spectrograms
	for(int chan=0; chan<16; chan++){
		for(int i=0; i<num_time_bins; i++){
			for(int j=0; j<num_freq_bins; j++){
				
				double temporary = spectrogram[chan]->GetBinContent(i,j); //get the current value of the bin
				if(number_of_entries->GetBinContent(i) >0.){temporary/=double(number_of_entries->GetBinContent(i));} //divide by the number of entries
				else temporary =1.; //unless there was nothing in that bin, in which case, we should just put a zero there
				spectrogram[chan]->SetBinContent(i,j,temporary); //set the bin content over again
				
				double temporary_background = spectrogram_background[chan]->GetBinContent(i,j); //get the current value of the bin
				if(number_of_entries_background->GetBinContent(i) >0.){temporary_background/=double(number_of_entries_background->GetBinContent(i));} //divide by the number of entries
				else temporary_background =1.; //unless there was nothing in that bin, in which case, we should just put a zero there
				spectrogram_background[chan]->SetBinContent(i,j,temporary_background); //set the bin content over again
				
				spectrogram_cleaned[chan]->SetBinContent(i,j,temporary-temporary_background); //subtract the two so you get a background subtracted answer
			}
		}
		//spectrogram_cleaned[chan] = 
		//spectrogram[chan]->Add(spectrogram_background,-1); //subtract the two
	}
	
	char title[150];
	
	TDatime plotstart(2011, 02, 15, 1, 55,0);
	int start_time = plotstart.Convert()-18000;
	TDatime plotstop(2011, 02, 15, 2, 02,0);
	int end_time = plotstop.Convert()-18000;
	
	double z_value_min = 1.E-2;
	double z_value_max = 3.E2;
	
	int plot=2;
	
	TCanvas *single = new TCanvas("","",1100,850);
	spectrogram[plot]->Draw("colz");
	spectrogram[plot]->GetYaxis()->SetRangeUser(60,400);
	spectrogram[plot]->GetXaxis()->SetRangeUser(start_time,end_time);
	spectrogram[plot]->GetZaxis()->SetRangeUser(z_value_min,z_value_max);
	spectrogram[plot]->GetYaxis()->SetTitle("Frequency (MHz)");
	spectrogram[plot]->GetXaxis()->SetTitle("UTC Time");
	spectrogram[plot]->GetZaxis()->SetTitle("Average Power Spectrum (dB(mV^{2}/MHz))");
	spectrogram[plot]->GetXaxis()->SetTitleOffset(1.1);
	spectrogram[plot]->GetYaxis()->SetTitleOffset(1.1);
	spectrogram[plot]->GetZaxis()->SetTitleOffset(1.1);
	spectrogram[plot]->GetXaxis()->SetTitleSize(0.045);
	spectrogram[plot]->GetYaxis()->SetTitleSize(0.045);
	spectrogram[plot]->GetZaxis()->SetTitleSize(0.045);
	spectrogram[plot]->GetXaxis()->SetLabelSize(0.045);
	spectrogram[plot]->GetYaxis()->SetLabelSize(0.045);
	spectrogram[plot]->GetZaxis()->SetLabelSize(0.045);
	gPad->SetLogz();
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.04);
	sprintf(title,"/n/home00/clark.2668/workspace/TrunkSolarFlares/results/today/%d.%d.%d_Paper_DynamicSpectrogram_%d_to_%d.pdf",year_now,month_now,day_now,start_time,end_time);
	single->SaveAs(title);
	delete single;
	
	TCanvas *single_cleaned = new TCanvas("","",1100,850);
	spectrogram_cleaned[plot]->Draw("colz");
	spectrogram_cleaned[plot]->GetYaxis()->SetRangeUser(60,400);
	spectrogram_cleaned[plot]->GetXaxis()->SetRangeUser(start_time,end_time);
	spectrogram_cleaned[plot]->GetZaxis()->SetRangeUser(z_value_min,z_value_max);
	spectrogram_cleaned[plot]->GetYaxis()->SetTitle("Frequency (MHz)");
	spectrogram_cleaned[plot]->GetXaxis()->SetTitle("UTC Time");
	spectrogram_cleaned[plot]->GetZaxis()->SetTitle("Average Power Spectrum (dB(mV^{2}/MHz))");
	spectrogram_cleaned[plot]->GetXaxis()->SetTitleOffset(1.1);
	spectrogram_cleaned[plot]->GetYaxis()->SetTitleOffset(1.1);
	spectrogram_cleaned[plot]->GetZaxis()->SetTitleOffset(1.1);
	spectrogram_cleaned[plot]->GetXaxis()->SetTitleSize(0.045);
	spectrogram_cleaned[plot]->GetYaxis()->SetTitleSize(0.045);
	spectrogram_cleaned[plot]->GetZaxis()->SetTitleSize(0.045);
	spectrogram_cleaned[plot]->GetXaxis()->SetLabelSize(0.045);
	spectrogram_cleaned[plot]->GetYaxis()->SetLabelSize(0.045);
	spectrogram_cleaned[plot]->GetZaxis()->SetLabelSize(0.045);
	gPad->SetLogz();
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.04);
	sprintf(title,"/n/home00/clark.2668/workspace/TrunkSolarFlares/results/today/%d.%d.%d_Paper_DynamicSpectrogramCleaned_%d_to_%d.pdf",year_now,month_now,day_now,start_time,end_time);
	single_cleaned->SaveAs(title);
	delete single_cleaned;
	
	TDatime plotstart_background(2011, 02, 11, 1, 55,0);
	int start_time_background = plotstart_background.Convert()-18000;
	TDatime plotstop_background(2011, 02, 11, 2, 02,0);
	int end_time_background = plotstop_background.Convert()-18000;
	
	TCanvas *single_background = new TCanvas("","",1100,850);
	spectrogram_background[plot]->Draw("colz");
	spectrogram_background[plot]->GetYaxis()->SetRangeUser(60,400);
	spectrogram_background[plot]->GetXaxis()->SetRangeUser(start_time_background,end_time_background);
	spectrogram_background[plot]->GetZaxis()->SetRangeUser(z_value_min,z_value_max);
	spectrogram_background[plot]->GetYaxis()->SetTitle("Frequency (MHz)");
	spectrogram_background[plot]->GetXaxis()->SetTitle("UTC Time");
	spectrogram_background[plot]->GetZaxis()->SetTitle("Average Power Spectrum (dB(mV^{2}/MHz))");
	spectrogram_background[plot]->GetXaxis()->SetTitleOffset(1.1);
	spectrogram_background[plot]->GetYaxis()->SetTitleOffset(1.1);
	spectrogram_background[plot]->GetZaxis()->SetTitleOffset(1.1);
	spectrogram_background[plot]->GetXaxis()->SetTitleSize(0.045);
	spectrogram_background[plot]->GetYaxis()->SetTitleSize(0.045);
	spectrogram_background[plot]->GetZaxis()->SetTitleSize(0.045);
	spectrogram_background[plot]->GetXaxis()->SetLabelSize(0.045);
	spectrogram_background[plot]->GetYaxis()->SetLabelSize(0.045);
	spectrogram_background[plot]->GetZaxis()->SetLabelSize(0.045);
	gPad->SetLogz();
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.04);
	sprintf(title,"/n/home00/clark.2668/workspace/TrunkSolarFlares/results/today/%d.%d.%d_Paper_DynamicSpectrogramBackground_%d_to_%d.pdf",year_now,month_now,day_now,start_time,end_time);
	single_background->SaveAs(title);
	delete single_background;
	
	TDatime t1i(2011, 02, 15, 1, 50,0);
	int t1i_int = t1i.Convert()-18000;
	TDatime t1f(2011, 02, 15, 2, 10,0);
	int t1f_int = t1f.Convert()-18000;
	
	TDatime t2i(2011, 02, 15, 2, 10,0);
	int t2i_int = t2i.Convert()-18000;
	TDatime t2f(2011, 02, 15, 2, 30,0);
	int t2f_int = t2f.Convert()-18000;
	
	TDatime t3i(2011, 02, 15, 2, 30,0);
	int t3i_int = t3i.Convert()-18000;
	TDatime t3f(2011, 02, 15, 2, 50,0);
	int t3f_int = t3f.Convert()-18000;
	
	int start_times[3] = {t1i_int,t2i_int,t3i_int};
	int end_times[3] = {t1f_int,t2f_int,t3f_int};
	
	TH2D* subplot[3];
	for(int i=0; i<3; i++){
		subplot[i] = new TH2D("","",num_time_bins,start_time_bin,stop_time_bin,num_freq_bins,start_freq_bin,stop_freq_bin); //the spectrogram   
		subplot[i] = (TH2D*) spectrogram[plot]->Clone();
		subplot[i]->GetXaxis()->SetTimeDisplay(1);
		subplot[i]->GetXaxis()->SetTimeFormat("%H:%M");
	}
	
	TCanvas *grid = new TCanvas("","",1100,1100);
	grid->Divide(1,3);
	//int canvas[9]={1,2,3,4,5,6,7,8,9};
	for(int i=0; i<3; i++){
		int start = start_times[i];
		int end = end_times[i];
		grid->cd(i+1);
		subplot[i]->Draw("colz");
		subplot[i]->GetYaxis()->SetRangeUser(60,400);
		subplot[i]->GetXaxis()->SetRangeUser(start,end);
		subplot[i]->GetZaxis()->SetRangeUser(z_value_min,z_value_max);
		subplot[i]->GetYaxis()->SetTitle("Frequency (MHz)");
		subplot[i]->GetXaxis()->SetTitle("UTC Time");
		subplot[i]->GetZaxis()->SetTitle("Average Power Spectrum (dB(mV^{2}/MHz))");
		subplot[i]->GetXaxis()->SetLabelSize(2*0.045);
		subplot[i]->GetYaxis()->SetLabelSize(2*0.045);
		subplot[i]->GetZaxis()->SetLabelSize(2*0.045);
		if(i==0){
			subplot[i]->GetXaxis()->SetTitleOffset(1);
			subplot[i]->GetYaxis()->SetTitleOffset(1);
			subplot[i]->GetZaxis()->SetTitleOffset(1);
			subplot[i]->GetXaxis()->SetTitleSize(2*0.045);
			subplot[i]->GetYaxis()->SetTitleSize(2*0.045);
			subplot[i]->GetZaxis()->SetTitleSize(2*0.045);
		}
		gPad->SetLogz();
		gPad->SetRightMargin(0.15);
		gPad->SetTopMargin(0.04);
		gPad->Update();
	}
	sprintf(title,"/n/home00/clark.2668/workspace/TrunkSolarFlares/results/today/%d.%d.%d_Paper_DynamicSpectrogram_ThreePanels.pdf)",year_now,month_now,day_now);
	grid->SaveAs(title);
	delete grid;
	
	for(int chan=0; chan<16; chan++) {delete spectrogram[chan]; delete spectrogram_background[chan]; delete spectrogram_cleaned[chan];}
	for(int i=0; i<3; i++) delete subplot[i];

}//close the main program
