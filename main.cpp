/*
 * Copyright (c) 2016, Gabriele Facciolo <gfacciol@gmail.com>
 *                     Nicola Pierazzo <nicolapierazzo@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>

#include "DCTdenoising.h"
#include "utils.hpp"
#include "demoutils.hpp"

using imgutils::pick_option;
using imgutils::read_image;
using imgutils::save_image;
using imgutils::Image;
using std::cerr;
using std::endl;
using std::move;
using std::vector;



int read_noise_lut(const char *name, vector<Image> &var_lut, vector<float> &bins_lut)
{
	FILE *f = fopen(name,"rt");

	// error checkings
	{
//		if (f == NULL) CV_Error( CV_StsObjectNotFound, "file not found");
//		if (mat.dims > 2) CV_Error( CV_StsBadSize, "function only handles 2D matrices" );
//		if (mat.channels() > 1) CV_Error ( CV_StsBadSize, "function only handles scalar matrices" );
	}

	int sz[2] = {0,0};
	// determine size
	{
		// number of rows
		int c;
		while (EOF != (c = fgetc(f))) if(c == '\n') ++sz[0];

		// number of cols
		rewind(f);
		for (sz[1] = 1; sz[1] < 2000; ++sz[1])
		{
			// we expect a floating point number
			float data;
			int success = fscanf(f,"%g",&data);

			// read next character
			c = fgetc(f);
			if ((c == '\n') || (c == EOF)) break;
		}
		rewind(f);
	}


	int nbins = sz[0];
	int dctsz = sqrt(sz[1] - 1);

	for (int b = 0; b < nbins; ++b)
	{
		float data;
		int success = fscanf(f, "%g", &data);
		bins_lut.push_back(data);

		Image bin_dct_var(dctsz,dctsz,1);

		for (int c = 0; c < dctsz; ++c)
		for (int r = 0; r < dctsz; ++r)
		{
			float data;
			int success = fscanf(f, "%g", &data);
			bin_dct_var.val(c,r) = data;
		}
		
		var_lut.push_back(move(bin_dct_var));

	}

  	fclose(f);
	return 1;
}


/**
 * @file   main.cpp
 * @brief  Main executable file. 
 *
 * @author Gabriele Facciolo
 * @author Nicola Pierazzo
 */
int main(int argc, char **argv) {
  // read DCT denoising options
  const bool  usage = static_cast<bool>(pick_option(&argc, argv, "h", nullptr));
  const int   dct_sz = atoi(pick_option(&argc, argv, "w", "8"));
  const char  *second_step_guide = pick_option(&argc, argv, "2", "");
  const char  *sigma_lut = pick_option(&argc, argv, "sigma_lut", "");
  const bool  no_second_step = static_cast<bool>(pick_option(&argc, argv, "1", NULL));
  const bool  no_first_step = second_step_guide[0] != '\0';
  const bool  adaptive_aggregation 
    = ! static_cast<bool>(pick_option(&argc, argv, "no_adaptive_aggregation", NULL));
  // read multiscaler options
  const char  *out_single = pick_option(&argc, argv, "single", "");
  const int   scales = atoi(pick_option(&argc, argv, "n", "4"));
  const float recompose_factor
    = static_cast<float>(atof(pick_option(&argc, argv, "c", ".5")));

  //! Check if there is the right call for the algorithm
  if (usage || argc < 2) {
    cerr << "usage: " << argv[0] << " sigma [input [output]] [-1 | -2 guide] "
         << "[-w patch_size (default 8)] [-c factor(.5)] [-n scales(4)] "
         << "[-sigma_lut sigma_lut_file] "
         << "[-single output_singlescale] [-no_adaptive_aggregation]" << endl;
    return usage ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  if (no_second_step && no_first_step) {
    cerr << "You can't use -1 and -2 together." << endl;
    return EXIT_FAILURE;
  }

#ifndef _OPENMP
  cerr << "Warning: OpenMP not available. The algorithm will run in a single" <<
       " thread." << endl;
#endif

  // read input
  Image noisy = read_image(argc > 2 ? argv[2] : "-");

  const float sigma = static_cast<float>(atof(argv[1]));
  vector<Image> sigma_dct_lut;
  vector<float> sigma_dct_lut_bins;
  if (sigma_lut[0] != '\0')
  {
	  printf("a\n");
	  printf("\%s\n", sigma_lut);
	  read_noise_lut(sigma_lut, sigma_dct_lut, sigma_dct_lut_bins);

//	  for (int i = 0; i < sigma_dct_lut_bins.size(); i++)
//	  {
//		  printf("%d %f\n", i, sigma_dct_lut_bins[i]);
//		  int d = sigma_dct_lut[0].rows();
//		  for (int r = 0; r < d; r++)
//		  {
//			  for (int c = 0; c < d; c++)
//				  printf("%f ", sigma_dct_lut[i].val(c,r));
//			  printf("\n");
//		  }
//	  }


//   printf("%dx%dx%d\n", sigma_dct_lut.rows(), sigma_dct_lut.columns(), sigma_dct_lut.channels());
	  printf("\%s\n", sigma_lut);
  }


  // generate the DCT pyramid 
  vector<Image> noisy_p = decompose(noisy, scales);
  vector<Image> guide_p, denoised_p;
  if (no_first_step) {
    Image guide = read_image(second_step_guide);
    guide_p = decompose(guide, scales);
  }

  // apply DCT denoising at each scale of the pyramid
  for (int layer = 0; layer < scales; ++layer) {
    // noise at the current scale is proportional to the number of pixels
    float s = sigma * sqrt(static_cast<float>(noisy_p[layer].pixels()) / noisy.pixels());
    if (!no_first_step) {
      Image guide = DCTdenoising(noisy_p[layer], s, dct_sz, sigma_dct_lut_bins, 
				                     sigma_dct_lut, adaptive_aggregation);
      guide_p.push_back(move(guide));
    }
    if (!no_second_step) {
      Image result =
          DCTdenoisingGuided(noisy_p[layer], guide_p[layer], s, dct_sz, sigma_dct_lut_bins,
					              sigma_dct_lut, adaptive_aggregation);
      denoised_p.push_back(move(result));
    } else {
      denoised_p.push_back(move(guide_p[layer]));
    }
  }

  // recompose pyramid
  if (strlen(out_single)) save_image(denoised_p[0], out_single);
  Image result = recompose(denoised_p, recompose_factor);

  save_image(result, argc > 3 ? argv[3] : "TIFF:-");

  return EXIT_SUCCESS;
}
