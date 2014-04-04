/*
Dipolar field function(s).

Copyright (c) 2014, Ilja Honkonen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the names of the copyright holders nor the names of their contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef BACKGROUND_B_DIPOLE_HPP
#define BACKGROUND_B_DIPOLE_HPP

#include "limits"

namespace background_B {

/*!
Returns magnetic field from a point dipole.

Returned value is a point value at given field position
from given dipole moment and position.
Assumes Vector_T is a 3d Eigen vector or similar.

Example usage:
\verbatim
const Eigen::Vector3d
	field = get_dipole_field(
		get_earth_dipole_moment<Eigen::Vector3d>(),
		{0, 0, 0},
		{6.3712e7, 0, 0} // at 10 Earth radii in +x
	);
\endverbatim
*/
template <class Vector_T> Vector_T get_dipole_field_0d(
	const Vector_T& dipole_moment,
	const Vector_T& dipole_position,
	const Vector_T& field_position
) {
	static_assert(
		Vector_T::SizeAtCompileTime == 3,
		"ERROR: Only 3 component vectors supported"
	);

	constexpr double permeability = 1.256637061e-6;

	const Vector_T r = field_position - dipole_position;

	bool far_enough = false;
	for (size_t i = 0; i < r.size(); i++) {
		if (fabs(r[i]) > std::numeric_limits<typename Vector_T::Scalar>::epsilon()) {
			far_enough = true;
			break;
		}
	}
	if (not far_enough) {
		return 2.0 / 3.0 * permeability * dipole_moment;
	}

	const double r3 = r.norm() * r.squaredNorm();
	const Vector_T
		r_unit = r / r.norm(),
		projected_dipole_moment = dipole_moment.dot(r_unit) * r_unit;
	return 1e-7 * (3 * projected_dipole_moment - dipole_moment) / r3;
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
