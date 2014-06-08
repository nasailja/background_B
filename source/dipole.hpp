/*
Dipolar field function(s).

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

#ifndef BACKGROUND_B_DIPOLE_HPP
#define BACKGROUND_B_DIPOLE_HPP

#include "limits"

#include "cubature.h"

namespace background_B {

/*!
Returns magnetic field from a point dipole.

Returned value is a point value at given field position
from given dipole moment and position.
Assumes Vector_T is a 3d Eigen vector or similar.

Example usage:
\verbatim
const Eigen::Vector3d
	field = get_dipole_field_0d(
		get_earth_dipole_moment<Eigen::Vector3d>(),
		{0, 0, 0},
		{6.3712e7, 0, 0} // at 10 Earth radii in +x
	);
\endverbatim
*/
template <class Vector_T> Vector_T get_dipole_field_0d(
	const std::vector<std::pair<Vector_T, Vector_T>>& dipole_moments_positions,
	const Vector_T& field_position
) {
	static_assert(
		Vector_T::SizeAtCompileTime == 3,
		"ERROR: Only 3 component vectors supported"
	);

	Vector_T ret_val;
	ret_val.setZero();

	constexpr double permeability = 1.256637061e-6;

	for (const auto dip_mom_pos: dipole_moments_positions) {

		const Vector_T r = field_position - dip_mom_pos.second;

		bool far_enough = false;
		for (size_t i = 0; i < r.size(); i++) {
			if (fabs(r[i]) > std::numeric_limits<typename Vector_T::Scalar>::epsilon()) {
				far_enough = true;
				break;
			}
		}
		if (not far_enough) {
			ret_val += 2.0 / 3.0 * permeability * dip_mom_pos.first;
			continue;
		}

		const double r3 = r.norm() * r.squaredNorm();
		const Vector_T
			r_unit = r / r.norm(),
			projected_dip_mom = dip_mom_pos.first.dot(r_unit) * r_unit;

		ret_val += 1e-7 * (3 * projected_dip_mom - dip_mom_pos.first) / r3;

	}

	return ret_val;
}


/*!
Returns magnetic field from a point dipole.

Returned value is the average value on the given line
from given dipole moment and position.
line_dimension is the dimension along which the line is.
Assumes Vector_T is a 3d Eigen vector or similar.

See http://ab-initio.mit.edu/wiki/index.php/Cubature
for help on maxEval, reqAbsError, reqRelError

Example usage:
\verbatim
TODO: const Eigen::Vector3d
	field = get_dipole_field(
		get_earth_dipole_moment<Eigen::Vector3d>(),
		{0, 0, 0},
		{6.3712e7, 0, 0} // at 10 Earth radii in +x
	);
\endverbatim
*/
template <class Vector_T> Vector_T get_dipole_field_1d(
	// pcubature requires non-const data
	std::vector<std::pair<Vector_T, Vector_T>>& dipole_moments_positions,
	const Vector_T& line_start,
	const typename Vector_T::Scalar line_length,
	const size_t line_dimension,
	const size_t maxEval = 0,
	const double reqAbsError = 0,
	const double reqRelError = 1e-6
) {
	static_assert(
		Vector_T::SizeAtCompileTime == 3,
		"ERROR: Only 3 component vectors supported"
	);

	if (line_dimension > 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Line dimension must be 0, 1 or 2: " << line_dimension
			<< std::endl;
		abort();
	}

	Vector_T line_end(line_start);
	line_end[line_dimension] += line_length;

	Vector_T integral, error;
	if (
		pcubature(
			Vector_T::SizeAtCompileTime,
			[](
				unsigned,
				const double* const x,
				void* external_data,
				unsigned,
				double* fval
			){
				const Eigen::Map<const Vector_T> pos(x);
				Eigen::Map<Vector_T> point_value(fval);

				point_value = get_dipole_field_0d<Vector_T>(
					*(static_cast<
						std::vector<std::pair<Vector_T, Vector_T>>*
					>(external_data)),
					pos
				);
				return 0;
			},
			static_cast<void*>(&dipole_moments_positions),
			Vector_T::SizeAtCompileTime,
			line_start.data(),
			line_end.data(),
			maxEval,
			reqAbsError,
			reqRelError,
			ERROR_LINF,
			integral.data(),
			error.data()
		) != 0
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Could not calculate line integral of a dipole from "
			<< line_start << " to " << line_end
			<< std::endl;
		abort();
	}

	std::cout << "returning" << std::endl << integral << std::endl;
	return integral;
}


/*!
Returns the Earth's dipole moment pointing in -Z direction.

Magnitude of the moment is 7.94e22 (A m^2), from
seismo.berkeley.edu/~rallen/eps122/lectures/L05.pdf
*/
template <class Vector_T> Vector_T get_earth_dipole_moment()
{
	return {0, 0, -7.94e22};	
}

/*!
Returns the Jupiter's dipole moment pointing in -Z direction.

Magnitude of the moment is 1.43e27 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_jupiter_dipole_moment()
{
	return {0, 0, -1.43e27};	
}

/*!
Returns the Saturn's dipole moment pointing in -Z direction.

Magnitude of the moment is 4.61e25 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_saturn_dipole_moment()
{
	return {0, 0, -4.61e25};	
}

/*!
Returns the Uranus' dipole moment pointing in -Z direction.

Magnitude of the moment is 3.97e24 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_uranus_dipole_moment()
{
	return {0, 0, -3.97e24};	
}

/*!
Returns the Neptune's dipole moment pointing in -Z direction.

Magnitude of the moment is 1.91e24 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_neptune_dipole_moment()
{
	return {0, 0, -1.91e24};	
}

} // namespace

#endif
