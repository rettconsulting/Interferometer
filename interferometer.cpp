//
// interferometer.cpp : virtual phase shifting interferometer with phase unwrapping
// Alex Martin - alex@rettc.com
//

#include <iostream>
#include <vector>
using namespace std;
#include "iniparser.h"
#include "cimg.h"
using namespace cimg_library;

int main()
{

	//configure from ini
	dictionary *profile;
	profile = iniparser_load("interferometer.ini");
	double lambda_start = iniparser_getdouble(profile, "Interferogram:Start_Wavelength", -1);
	double lambda_step = iniparser_getdouble(profile, "Interferogram:Wavelength_Step", -1);
	int number_of_steps = iniparser_getint(profile, "Interferogram:Number_Steps", -1);
	int mu = iniparser_getint(profile, "Interferogram:Resolution_Padding", -1);
	double z_part_thickness = iniparser_getdouble(profile, "Interferogram:Part_Thickness", -1);
	iniparser_freedict(profile);

	double C = 299792500; // speed of light, meters per second
   double Pi = 3.14159; // apple pie

	// convert nanometers to meters
	lambda_start /= 1000000000;
	lambda_step /= 1000000000;

	cout << "Start Lambda = " << lambda_start << " meters" << endl;
	cout << "Start Frequency = " << C/lambda_start << " hz" << endl;

	cout << "Lambda Step = " << lambda_step << " meters" << endl;
	cout << "Frequency Step = " << (C/lambda_start)/lambda_start*lambda_step << " hz" << endl;
	double frequency_step = ((C/lambda_start)/lambda_start)*lambda_step;

	//add 0.5 so that cast to int rounds properly
//	int number_samples = (int)(((lambda_end-lambda_start)/lambda_step)+0.5);
	
	int number_samples = number_of_steps;

	cout << "Number of Steps = " << number_samples << endl << endl;

	// setup constants
	const unsigned char grey[] = {128,128,128};
	const unsigned char white[] = { 255,255,255 }; 
	const unsigned char red[] = {255,0,0};
	const unsigned char blue[] = {0,0,255};
	const unsigned char black[] = {0,0,0};

	double lambda = lambda_start;
	double frequency = C/lambda_start;

	// generate interferogram
	// init image
	
	int height = 128;
	int width = 128;

	//initialize surfaces, initialize to all zeros
	CImg<double> surface_transmission(width, height,1, 1,0);
	CImg<double> surface_part_input(width, height,1, 1,0);
	CImg<double> surface_part_output(width, height,1, 1,0);
	CImg<double> surface_reflection(width, height,1, 1,0);

	//initialize the image for the interferogram
	CImg<double> interferogram(width, height,1, 1,0);
	CImg<double> ref_interferogram(width, height,1, 1,0);

// #########################################################
// #   surface parameters: generate 4 surface geometry:    #
// #########################################################

// 1) flat transmission
// 2) part side one
// 3) part side two
// 4) flat reflection

   //init surface generation parameter variables.

	double k_sphere = 0; 
	int x_offset = 0;
	int y_offset = 0;

// generate surface_transmission

	//this surface is already initialized to zero
	//this is our z axis reference z = 0
   for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			//surface_transmission(x,y) += ( sin((double)20*y/height) + cos((double)30*x/width) )/(2*Pi*10000000);
			//surface_transmission(x,y) += ( y/(6*Pi*100000000));
			//surface is perfectly flat.
      }
	}

// generate surface_part_input (closest to zero)
	//space as 1st order [(1+mu) * z_part_thickness]
	surface_part_input += z_part_thickness*(1+mu);
	surface_part_input /= 1000; //convert to meters 
	// power
	k_sphere = -1000; 
	x_offset = width/2;
	y_offset = height/2;

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			//add power
			surface_part_input(x,y) += k_sphere*( (pow((double)(x-x_offset)/(double)width , 2) + pow((double)(y-y_offset)/(double)height,2))/1000000000 );
		}
	}

// generate surface_part_output
	//space as 1st order + 0th order = [(1+mu) * z_part_thickness] + [z_part_thickness]
	surface_part_output += z_part_thickness*(1+mu) + z_part_thickness;
	surface_part_output /= 1000; //convert to meters
	// power
	k_sphere = 2500; 
	x_offset = width/3; // offset from center to left
	y_offset = height/2;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			//add power
			//surface_part_output(x,y) += k_sphere*( (pow((double)(x-x_offset)/(double)width , 2) + pow((double)(y-y_offset)/(double)height,2))/1000000000 );
			surface_part_output(x,y) += (sin((double)5*y/height) + sin((double)3*x/width))/(2*Pi*100000);
		}
	}


//generate surface_reflection
	//space as 1st order + 0th order + 2nd order = 
	//[(1*mu)*z_part_thickness] + [z_part_thickness] + [(1*mu)*(1*mu)*z_part_thickness]
	surface_reflection += z_part_thickness*(1+mu) + z_part_thickness + z_part_thickness*(1+mu)*(1+mu);
	surface_reflection /= 1000; // convert to meters 
	// x tilt
	double x_slope = 30;
	// y_tilt
	double y_slope = 7;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			//add tilt
			//surface_reflection(x,y) += (y_slope*y + x_slope*x)/1000000000;
			surface_reflection(x,y) += ( sin((double)20*y/height) + cos((double)30*x/width) )/(2*Pi*10000000);
		}
	}

// output approximate spacing (spacing before adding figure)
	cout << "Transmission to S1 Length = " << (z_part_thickness*(1+mu))/1000 << " meters" << endl;
	cout << "S1 to S2 Length = " << (z_part_thickness)/1000 << " meters" << endl;
	cout << "S2 to Reflection = " << (z_part_thickness*(1+mu)*(1+mu))/1000 << " meters" <<  endl;
	cout << "Transmission to Reflection = "
		<< (z_part_thickness*(1+mu)*(1+mu))/1000 
		+ (z_part_thickness*(1+mu))/1000 
		+ (z_part_thickness)/1000 << " meters" <<  endl;

	//make temp images for display
	CImg<double> surface_transmission_temp = surface_transmission;
	CImg<double> surface_part_input_temp = surface_part_input;
	CImg<double> surface_part_output_temp = surface_part_output;
	CImg<double> surface_reflection_temp = surface_reflection;

	//display 4 surfaces, flat, s1, s2, tilted reflection.
	CImgList<double> surfaces(surface_transmission_temp.normalize(0,255), surface_part_input_temp.normalize(0,255), surface_part_output_temp.normalize(0,255), surface_reflection_temp.normalize(0,255));
	CImgDisplay surfaces_disp(surfaces, "Initial Surface Display");
	surfaces.display(surfaces_disp);
	
	// init inteferogram display
	CImgDisplay interferogram_disp; 

	// initialize interferogram image storage lists
	CImgList<double> dataset;
	CImgList<double> ref_dataset;

	// seed for rand() for lambda jitter
	srand ( (int)time(NULL) );
	
	for (int i = 0; i < number_samples; i++) {

		CImgList<double> interferograms(interferogram.normalize(0,1),ref_interferogram.normalize(0,1));
		interferogram_disp.set_title("Interferogram Display");
		interferograms.display(interferogram_disp);
      
		// iterate wavelength
		lambda += lambda_step;
		frequency += frequency_step;
//		cout << "frequency = " << frequency << " lambda = " << C/frequency <<endl;

		// ADD JITTER HERE
		double jitter = (((double)rand()/((double)RAND_MAX)-0.5)/10000000)+1;
      // screw jitter
		jitter = 1;
//		cout << "jitter = " << jitter-1 << endl;

		// generate new phase map
		for (int y = 0; y < interferogram.height(); y++) {
			for (int x = 0; x < interferogram.width(); x++) {
				// sum up cavity differences
				// b-a, c-a, d-a
				interferogram(x,y) = pow( cos((4*Pi*frequency/(C*jitter))*(surface_part_input(x,y) - surface_transmission(x,y)) ), 2 );
				interferogram(x,y) += pow( cos((4*Pi*frequency/(C*jitter))*(surface_part_output(x,y) - surface_transmission(x,y)) ), 2 );
				interferogram(x,y) += pow( cos((4*Pi*frequency/(C*jitter))*(surface_reflection(x,y) - surface_transmission(x,y)) ), 2 );
				// c-b, d-b
				interferogram(x,y) += pow( cos((4*Pi*frequency/(C*jitter))*(surface_part_output(x,y) - surface_part_input(x,y)) ), 2 );
				interferogram(x,y) += pow( cos((4*Pi*frequency/(C*jitter))*(surface_reflection(x,y) - surface_part_input(x,y)) ), 2 );
				// d-c
				interferogram(x,y) += pow( cos((4*Pi*frequency/(C*jitter))*(surface_reflection(x,y) - surface_part_output(x,y)) ), 2 );
				// sum up reference cavity
				ref_interferogram(x,y) = pow( cos((4*Pi*frequency/(C*jitter))*(surface_reflection(x,y) - surface_transmission(x,y)) ), 2 );
			}
		}
		// append interferogram to dataset
		dataset.push_back(interferogram);
  		ref_dataset.push_back(ref_interferogram);

    } // end data collection


//#######################################################
//#				GET FFTs FOR OPL SPECTRUM				#
//#######################################################

	cout << "\nGetting FFTs:";

	CImg<double> fft_set_re(width, height, dataset.size(), 1, 0);
	CImg<double> fft_set_im(width, height, dataset.size(), 1, 0);
	CImg<double> opl_mag_set(width, height, dataset.size()/2, 1, 0);

	// iterate through x and y coords
	for (int x = 0; x < width; x++) {
		cout << ".";
		for (int y = 0; y < height; y++) {
			// compute fft of one x,y point out of all data along z
			// 1d data set (dx, dy, dz, dc, init)
			CImg<double> opl_spectrum(dataset.size(),1,1, 1,0);
			// populate dataset over samples
			for (unsigned int z = 0; z < dataset.size(); z++) {
				opl_spectrum(z) = dataset[z](x,y);

			}

			// get fft (duh)
			CImgList<float> F = opl_spectrum.get_FFT(); 
			
			//generate OPL spectrum
			for (unsigned int z = 0; z < dataset.size(); z++) {
				fft_set_re(x,y,z) = F[0](z,0);
				fft_set_im(x,y,z) = F[1](z,0);
			}
			
			//
			CImg<double> opl_mag = ( (F[0].get_pow(2) + F[1].get_pow(2)).sqrt() +1 ).log(); 
			//erase dc term
			opl_mag.crop(0,0,0,0,opl_mag.width()/2,0,0,0);
			opl_mag(0,0) = 0;	
			// copy x,y into z
			for (int z = 0; z < opl_mag_set.depth(); z++) {
				opl_mag_set(x,y,z) = opl_mag(z,0);
			}
		}
	}

	cout << endl << "Done Getting FFTs!\n" << endl;



//#######################################################
//#         	CREATE OPL PEAKS DATASET	            #
//#######################################################

	cout << "Finding OPL Peaks:";

	// create datum struct
	struct datum {
		int index;
		double value;
	};

	//create struct for peaks at a given x,y coord
	struct x_y_peak_set {
		int x; 
		int y;
		vector<datum> peaks;
	};

	// create vector of x_y_peaks
	vector<x_y_peak_set> opl_mag_peak_set;

	// iterate through x,y
	for (int x = 0; x < opl_mag_set.width(); x++) {	
		cout << ".";
		for (int y = 0; y < opl_mag_set.width(); y++) {
	
			// create image from crop along z axis at x,y
			CImg<double> opl_mag( opl_mag_set.get_crop(x,y,0,0,x,y,opl_mag_set.depth()-1,0) );

			// initialize peak_set
			struct x_y_peak_set peak_set;
			// set x,y coords for peak_set			
			peak_set.x = x;
			peak_set.y = y;

			// create derivatives vector
			vector<datum> opl_derivatives;

			int d_step = 3; // how many points ahead to take derivative TODO: PUT INT INI
			double d_tolerance = 1; // tolerance TODO: PUT INT INI

			// take derivatives
			for (int i = 0; i< (opl_mag.depth() - d_step); i++) {
				struct datum derivative; 
				derivative.index = i;
				derivative.value = opl_mag(0,0,i)-opl_mag(0,0,i+d_step);
				opl_derivatives.push_back(derivative);
			}
			
			double mag_tolerance_percent = 0.75; // 75%
			double mag_min = 0;
			double mag_max = opl_mag.min_max(mag_min); // get maximum magnitude
			double mag_tolerance = mag_max * mag_tolerance_percent; // set peak tolerance level

			// tolerance out small peaks, store peaks
			for (unsigned int i = 0; i < opl_derivatives.size(); i++) {
				// check derivative tolerance, and magnitude tolerance
				if (abs(opl_derivatives[i].value) > d_tolerance && opl_mag(i,0) > mag_tolerance) {
					// check neighbors of actual magnitude to make sure you are on peak
					if (opl_mag(i-1,0) < opl_mag(i,0) && opl_mag(i+1,0) < opl_mag(i,0)) {
						struct datum peak;
						peak.index = i;
						peak.value = opl_mag(i,0);
						peak_set.peaks.push_back(peak);
					
					}// end of neighbor check				
				}// 
			} // end of final peak finder (i) loop
			// append peak set for given x,y coord to ugly storage doo da
			opl_mag_peak_set.push_back(peak_set);
		}
	}

// fudging some peak finding errors by taking closest integer index to average across x,y peaksets.

	// get average number of peaks

	// now we find averages and get best integer case

	double average_num_peaks = 0;
	int integer_num_peaks = 0;

	for (int y = 0; y < height; y++) {

		for (int x = 0; x < width; x++) {
			average_num_peaks += opl_mag_peak_set[y*width+x].peaks.size();
		}
	}
	average_num_peaks /= width*height;

	// round to int
	integer_num_peaks = (int)floor(average_num_peaks+0.5);
	// initialize peak set for index values
	vector<int> average_peak_set;

	// make sure peak finder output makes sense, 4 surface cavity = 6 peaks
	if (integer_num_peaks == 6) {

		// get average index at given peaks
		for (int z = 0; z < integer_num_peaks; z++) {
			// init 
			double averaged_peak_index = 0;
			int integer_peak_index = 0;
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					averaged_peak_index += opl_mag_peak_set[y*width+x].peaks[z].index;
				}
			}
			// normalize
			averaged_peak_index /= width*height;
	//		cout << "averaged_peak_index = " << averaged_peak_index << endl;
			// get integer index value
			integer_peak_index = (int)floor(averaged_peak_index+0.5);
	//		cout << "integer_peak_index = " << integer_peak_index << endl;
			average_peak_set.push_back(integer_peak_index);
		}	

		cout << endl << "Done Finding OPL Peaks." << endl;
		cout << "Number of Peaks: " << average_peak_set.size() << endl;
	} 
	// if peak finder failed, debug display, and setup for safe finish of execution
	else {
		cout << endl << "Cannot Find Peaks!" << endl;
		cout << "Number of Peaks: " << average_num_peaks << endl;

		int zoom = 2;
		// init image for graphing OPL spectrum values  
		CImg<unsigned char> debug_fft_mag(zoom*dataset.size()/2,zoom*dataset.size()/2,1,3,0);
		// init image display for graphing OPL spectrum values
		CImgDisplay debug_fft_disp(debug_fft_mag.resize_tripleXY(),"OPL Spectrum");
		// draw graph
		debug_fft_mag.draw_graph(opl_mag_set.get_crop(0,0,0,0,0,0,opl_mag_set.depth()-1,0),white,1,1,4,0,0);
		// display image
		debug_fft_mag.display(debug_fft_disp);
	}


//#######################################################
//#         		SURFACE PHASE EXTRACTION	        #
//#######################################################

	cout << endl << "Calculating Cavity Phase Maps:" << endl;

	// create image list for phase maps

	CImgList<double> phase_maps;

	for (unsigned int z = 0; z < average_peak_set.size(); z++) {
		// init temporary phase image
		CImg<double> opl_arg(width, height, 1,1,0);
		for (int y = 0; y < height; y++) {
			//cout << ".";
			for (int x = 0; x < width; x++) {
				double im = fft_set_im(x,y,average_peak_set[z]);
				double re = fft_set_re(x,y,average_peak_set[z]);
				// calculate argument of phase
				opl_arg(x,y) = atan2(im,re);
			}
		}
		// append temporary phase image to image list
		phase_maps.push_back(opl_arg);
		cout << "Phase Map " << z+1 << " calculated." << endl;
	}

	// display phase maps
	CImgDisplay phase_disp(phase_maps, "Phase Display");
	//cout << "phase_maps.size() = " << phase_maps.size() << endl;
	if (phase_maps.size() > 0) {
		phase_maps.display(phase_disp);
	}

	cout << "Cavity Phase Map Calculation Complete." << endl << endl;

//#######################################################
//#         	UNWRAPPED SURFACE EXTRACTION	        #
//#######################################################

	cout << endl << "Unwrapping Cavities:" << endl;

	// create image list for phase maps

	CImgList<double> unwrapped_cavities;

	// create artificial value as for bad data indicator
	double NIL_DATA = 1024; // large positive number
	double MAX_DELTA = 0.75*2*Pi; // 75% of 2Pi

	//TODO filter phase maps, add NIL_DATA per some algorithm
	//right now we are using perfect virtual data thus perfect phase maps

	//why is phase_maps.size() returning a double????
	//cout << "phase_maps.size() = " << phase_maps.size() << endl;
	for (unsigned int i = 0; i < (int)phase_maps.size(); i++) {

		// init temporary phase image to NIL_DATA value
		CImg<double> unwrapped_cavity(width, height, 1,1,NIL_DATA);

		// create copy of real phase map to work with
		CImg<double> wrapped_cavity(phase_maps[i]);

		// add NIL_DATA border to wrapped_cavity
		// add top and bottom row
		for (int x = 0; x < width; x++) {
			wrapped_cavity(x,0) = NIL_DATA;
			wrapped_cavity(x,height-1) = NIL_DATA;
		}
		// add left and right column
		for (int y = 0; y < width; y++) {
			wrapped_cavity(0,y) = NIL_DATA;
			wrapped_cavity(width-1,y) = NIL_DATA;
		}

//#################################################################

		int  y, x, start_x, x1, x2;
		bool done, leak_found;
		double last, offset, start_offset, delta;
		int num_pixels=0;

		// Start the algorithm at the passed coordinates with no offset
      // Starting at the center
		y = width/2;
		x = height/2;
		x1 = 0;
		x2 = 0;
		offset = 0.0;
		done = false;
		int * search_stack = new int[50000];
		int * sp;
		sp = search_stack;
		int count = 0;

		while (!done) {
			count++;
			// Proceed from new starting position
			start_x = x;
			start_offset = offset;
			last = wrapped_cavity(x,y);

			unwrapped_cavity(x,y) = wrapped_cavity(x,y) + offset;  
			++num_pixels;

			// Go left as far as possible connecting phases as we go
			--x;
			while (wrapped_cavity(x,y) != NIL_DATA) {
					delta = wrapped_cavity(x,y) - last;
					if      (delta > +MAX_DELTA) {
						offset -= 2*Pi;
					}
					else if (delta < -MAX_DELTA) {
						offset += 2*Pi;
					}
					unwrapped_cavity(x,y) = wrapped_cavity(x,y) + offset;  ++num_pixels;
					last = wrapped_cavity(x,y);
					--x;
			}
			x1 = x + 1;    // Record left end segment

			// Reset to starting position
			x = start_x;
			offset = start_offset;
			last = wrapped_cavity(x,y);

			// Go right as far as possible connecting phases as we go
			++x;
			while (wrapped_cavity(x,y) != NIL_DATA) {
				delta = wrapped_cavity(x,y) - last;
				if      (delta > +MAX_DELTA) offset -= 2*Pi;
				else if (delta < -MAX_DELTA) offset += 2*Pi;
				unwrapped_cavity(x,y) = wrapped_cavity(x,y) + offset;  ++num_pixels;
				last = wrapped_cavity(x,y);
				++x;
			}
			x2 = x - 1;  // Record right end of segment
			
			// ******************************************************************
			// Search above and below for an available pixel to move to
			leak_found = false;
			while (!leak_found) {
				// ************* Look Above ************
				// Search pixels above current for hole going from left to right
				x = x1;
				while ( ( (wrapped_cavity(x,y+1) == NIL_DATA)       // No hole above
						|| (unwrapped_cavity(x,y+1) != NIL_DATA) )  // Pixel already cnnctd
						&& (x <= x2) ) {                 // But not to right extreme
					++x;
				}

				if (x > x2) {
					// ************* Look Below ************
					// Search pixels below for hole going from left to right
					x = x1;
					while (((wrapped_cavity(x,y-1) >= NIL_DATA) || (unwrapped_cavity(x,y-1) <  NIL_DATA)) && (x <= x2) ) {
						// No hole below // But not to right extreme  // Pixel already cnnctd
						++x;
					}

					if (x > x2) {   // No empty pixel below this segment
						// If no more unsearched segments then all done
						if (sp == search_stack) { done = true; break; }
						// Restore most recent unsearched segment from the stack
						--sp; x2=*sp; --sp; x1=*sp; --sp; y=*sp;    // Pop
						// Re-establish phase offset
						offset = unwrapped_cavity(x,y) - wrapped_cavity(x,y);
					}
					// ************* Leak downward ************
					else {
						// Save the current segment (y,x1,x2) on stack
						*sp=y; ++sp; *sp=x1; ++sp; *sp=x2; ++sp;    // Push
						// Re-establish phase offset
						offset = unwrapped_cavity(x,y) - wrapped_cavity(x,y);
						// Check for phase discontinuity going down
						last = wrapped_cavity(x,y);
						--y;
						delta = wrapped_cavity(x,y) - last;
						if      (delta > +MAX_DELTA) offset -= 2*Pi;
						else if (delta < -MAX_DELTA) offset += 2*Pi;
						// Exit leak search loop
						leak_found = true;
					}
				}
				// ************* Leak upward ************
				else {
					// Save the current segment (y,x1,x2) on stack
					*sp=y; ++sp; *sp=x1; ++sp; *sp=x2; ++sp;       // Push
					// Re-establish phase offset
					offset = unwrapped_cavity(x,y) - wrapped_cavity(x,y);
					// Check for discontinuity going up
					last = wrapped_cavity(x,y);
					++y;
					delta = wrapped_cavity(x,y) - last;
					if      (delta > +MAX_DELTA) offset -= 2*Pi;
					else if (delta < -MAX_DELTA) offset += 2*Pi;
					// Exit leak search loop
					leak_found = true;
				}
				if (count > 10000) {
					break;
				}
			}
		}

//#################################################################

		delete [] search_stack;

		// replace NIL_DATA with top and bottom row of zeros
		for (int x = 0; x < width; x++) {
			unwrapped_cavity(x,0) = 0.0;
			unwrapped_cavity(x,height-1) = 0.0;
		}
		// replace NIL_DATA with left and right column of zeroes
		for (int y = 0; y < width; y++) {
			unwrapped_cavity(0,y) = 0.0;
			unwrapped_cavity(width-1,y) = 0.0;
		}

		//TEMP normalize surfaces
//		unwrapped_cavity.normalize(0,255);
		// push onto stack
		unwrapped_cavities.push_back(unwrapped_cavity);
		cout << "Unwrapped Cavity Map " << i+1 << " calculated." << endl;
	}
	//cout << "unwrapped_cavities.size() = " << unwrapped_cavities.size() << endl;

	// copy data for surface calculations below
	
	CImg<double> BC,AB,AC,CD,BD,AD;
	
	BC = unwrapped_cavities[0];
	AB = unwrapped_cavities[1];
	AC = unwrapped_cavities[2];
	CD = unwrapped_cavities[3];
	BD = unwrapped_cavities[4];
	AD = unwrapped_cavities[5];


	// display unwrapped cavity maps
	
	unwrapped_cavities[0].normalize(0,255);
	unwrapped_cavities[1].normalize(0,255);
	unwrapped_cavities[2].normalize(0,255);
	unwrapped_cavities[3].normalize(0,255);
	unwrapped_cavities[4].normalize(0,255);
	unwrapped_cavities[5].normalize(0,255);
	
	CImgDisplay unwrapped_cavities_disp(unwrapped_cavities, "Unwrapped Cavity Display");
	if (unwrapped_cavities.size() > 0) {
		//cout << "(unwrapped_cavities.size() > 0)" << endl;
		unwrapped_cavities.display(unwrapped_cavities_disp);
	}

// subtract cavities to get original surfaces
	
	// cavities: A = transmission flat, B = part input, C = part output, D = reflection flat
	
	// unwrapped cavities are ordered by length, shortest first
	
	// unwrapped cavity 1 = BC => part cavity
	// unwrapped cavity 2 = AB => trans to part input cavity
	// unwrapped cavity 3 = AC => trans to part output cavity
	// unwrapped cavity 4 = CD => part output to reflection cavity
	// unwrapped cavity 5 = BD => part input to reflection cavity
	// unwrapped cavity 6 = AD => trans to reflection cavity
	
	// init surface image 
	CImg<double> unwrapped_surface(width, height, 1,1,0);
	// init surface list
	CImgList<double> unwrapped_surfaces;

	// solve for transmission surface	
	// this is already set to zero from out reference transmission flat;
	unwrapped_surface.normalize(0,255);
	unwrapped_surfaces.push_back(unwrapped_surface);
	// solve for part input surface
	unwrapped_surface = AB; // since A is reference transmission flat and is known
	unwrapped_surface.normalize(0,255);
	unwrapped_surfaces.push_back(unwrapped_surface);
	// solve for part output surface
	unwrapped_surface = -CD + AD; // again where A is zero
	unwrapped_surface.normalize(0,255);
	unwrapped_surfaces.push_back(unwrapped_surface);
	// solve for reflection surface
	unwrapped_surface = -AD;
	unwrapped_surface.normalize(0,255);
	unwrapped_surfaces.push_back(unwrapped_surface);


	// display unwrapped surface maps6
	CImgDisplay unwrapped_surfaces_disp(unwrapped_surfaces, "Unwrapped Surface Display");
	if (unwrapped_surfaces.size() > 0) {
		//cout << "(unwrapped_surfaces.size() > 0)" << endl;
		unwrapped_surfaces.display(unwrapped_surfaces_disp);
	}

	cout << endl << "Surface Phase Unwrapping Complete." << endl << endl;

//#######################################################
//#         		calculate surface distances	        #
//#######################################################

   cout << "Optical Path Lengths of Fundamental Cavities: " << endl;

   // based on fft spacing, absolute distances of surfaces are given per OPL spectrum
   // convert from wavelength to frequency
   double neu_start = C/(lambda_start);
   // convert from wavelength to frequency step
   double neu_step = neu_start*lambda_step/lambda_start;
   //iterate through average OPL peak set and print OPL lengths in meters
   cout << "OPL resolution = " << C/(number_of_steps*neu_step*4) << " meters." << endl;

   for (unsigned int i = 0; i < (unsigned int)average_peak_set.size(); i++) {
      // print optical path length
      double OPL = (average_peak_set[i]*C) / (number_of_steps*neu_step*4);
      cout << "Peak " << i << " = " << OPL << " meters." << endl;
   }
   cout << endl << "Compare to original Path Lengths:" << endl;

   // output approximate spacing (spacing before adding figure)
   cout << "Transmission to S1 Length = " << (z_part_thickness* (1 + mu)) / 1000 << " meters" << endl;
   cout << "S1 to S2 Length = " << (z_part_thickness) / 1000 << " meters" << endl;
   cout << "S2 to Reflection = " << (z_part_thickness * (1 + mu) * (1 + mu)) / 1000 << " meters" << endl;
   cout << "Transmission to Reflection = "
	   << (z_part_thickness * (1 + mu) * (1 + mu)) / 1000
	   + (z_part_thickness * (1 + mu)) / 1000
	   + (z_part_thickness) / 1000 << " meters" << endl;

   cout << endl << "Fourier Transform Frequency Scanning Interferometry demonstration complete." << endl << endl;

   cout << "Based on work by Leslie Deck of the Zygo Corporation. See his paper:" << endl << endl;
   cout << "https://www.osapublishing.org/ao/abstract.cfm?uri=ao-42-13-2354" << endl;

   cout << endl << "Copyright © Alex Martin - alex@rettc.com - http://www.rettc.com" << endl << endl;


//#######################################################
//#		OPL PEAK / SPECTRA EXPLORER - LOOP TILL QUIT	#
//#######################################################

	int zoom = 2;
	// init image for graphing OPL spectrum values  
	CImg<unsigned char> fft_mag(zoom*dataset.size()/2,zoom*dataset.size()/2,1,3,0);
	// init image display for graphing OPL spectrum values
	CImgDisplay fft_disp(fft_mag,"OPL Spectrum");
	// draw graph
	fft_mag.draw_graph(opl_mag_set.get_crop(0,0,0,0,0,0,opl_mag_set.depth()-1,0),white,1,1,4,0,0);
	// display image
	fft_mag.display(fft_disp);

	// loop until interferogram display is closed
	while (!interferogram_disp.is_closed() && !surfaces_disp.is_closed() && !fft_disp.is_closed()) {
//	while (!interferogram_disp.is_closed() && !surfaces_disp.is_closed() && !fft_disp.is_closed() && !phase_disp.is_closed()) {
    interferogram_disp.wait();
      if (interferogram_disp.button() && interferogram_disp.mouse_y()>=0 && interferogram_disp.mouse_x()>=0 && interferogram_disp.mouse_x()< width) {
        const int y = interferogram_disp.mouse_y();
        const int x = interferogram_disp.mouse_x();
		//init display 
		fft_mag.fill(0);
		fft_mag.draw_graph(opl_mag_set.get_crop(x,y,0,0,x,y,opl_mag_set.depth()-1,0),white,1,1,4,0,0);
		int num_peaks = opl_mag_peak_set[y*width+x].peaks.size();
		for (int i = 0; i < num_peaks; i++) {
			fft_mag.draw_line(zoom*opl_mag_peak_set[y*width+x].peaks[i].index,0,zoom*opl_mag_peak_set[y*width+x].peaks[i].index, fft_mag.height(), blue);
			fft_mag.draw_text(10,20*i, "i = %i",red, NULL, 1, 20, opl_mag_peak_set[y*width+x].peaks[i].index);
			fft_mag.draw_text(80,20*i, "val = %4.5f",red, NULL, 1, 20, (float)opl_mag_peak_set[y*width+x].peaks[i].value);
		}
		fft_mag.draw_text(fft_mag.width()/2,fft_mag.height()/2, "(%i,%i)\n%i",red, NULL, 1, 36, x,y,num_peaks);
		fft_mag.display(fft_disp);
      }
    }
	return 0;
}
