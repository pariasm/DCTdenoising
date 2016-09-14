/*
 * Original code Copyright (c) 2010, Guoshen Yu <yu@cmap.polytechnique.fr>,
 *                                   Guillermo Sapiro <guille@umn.edu>
 * Modified code Copyright (c) 2015, Gabriele Facciolo <gfacciol@gmail.com>
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


/*--------------------------- DCTdenoising  -------------------------*/
// Detect corresponding points in two images with the ASIFT method.
// Copyright, Guoshen Yu, Guillermo Sapiro, 2010.
// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
/*---------------------------------------------------------------------------*/

#ifndef DCTDENOISING_DCTDENOISING_HPP
#define DCTDENOISING_DCTDENOISING_HPP

#include "Image.hpp"

imgutils::Image DCTdenoising(const imgutils::Image &noisy, float sigma,
                             int dct_size, int nthreads = 0);
imgutils::Image DCTdenoisingGuided(const imgutils::Image &noisy,
                                   const imgutils::Image &guide,
                                   float sigma, int dct_size, int nthreads = 0);

#endif  // DCTDENOISING_DCTDENOISING_HPP
