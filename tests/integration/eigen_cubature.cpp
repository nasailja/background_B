/*
Tests numerical integration using cubature with Eigen vectors.

Copyright 2014 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "array"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "limits"

#include "Eigen/Core"

#include "cubature.h"

int main()
{
	// integrate f(x) = x*y*y*z*z*z
	const Eigen::Vector3d
		start(-3, -2, -1),
		end(4, 3, 2);

	double integral, error;
	if (
		pcubature(
			1,
			[](
				unsigned,
				const double* const given_r,
				void*,
				unsigned,
				double* fval
			){
				const Eigen::Map<const Eigen::Vector3d> r(given_r);
				*fval = r[0] * r[1] * r[1] * r[2] * r[2] * r[2];
				return 0;
			},
			NULL,
			3,
			start.data(),
			end.data(),
			0,
			0,
			1e-3,
			ERROR_LINF,
			&integral,
			&error
		) != 0
	) {
		abort();
	}

	const double
		i = integral,
		e = 153.125;
	if (std::fabs(i - e) > 1e-9 * std::max(std::fabs(i), std::fabs(e))) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Too large relative error when integrating from\n" << start
			<< "\nto\n" << end << "\nexact: " << e << ", numerical: " << i
			<< ", error: " << error
			<< std::endl;
		abort();
	}


	// integrate f(x,y,z) = |r|, r.r, r.r|r|
	Eigen::Vector3d integral3d, error3d;
	if (
		pcubature(
			3,
			[](
				unsigned,
				const double* const given_r,
				void*,
				unsigned,
				double* fval
			){
				const Eigen::Map<const Eigen::Vector3d> r(given_r);
				Eigen::Map<Eigen::Vector3d> value(fval);
				value[0] = r.norm();
				value[1] = r.dot(r);
				value[2] = r.norm()*r.dot(r);
				return 0;
			},
			NULL,
			3,
			start.data(),
			end.data(),
			0,
			0,
			1e-9,
			ERROR_LINF,
			integral3d.data(),
			error3d.data()
		) != 0
	) {
		abort();
	}

	return EXIT_SUCCESS;
}
