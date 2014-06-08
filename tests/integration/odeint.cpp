/*
Tests numerical integration using boost odeint.

Copyright 2014 Ilja Honkonen
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

#include "array"
#include "boost/numeric/odeint.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"

using integral_type = std::array<double, 1>;

int main()
{
	/*
	Integrate functions with one argument from start to end
	*/

	const auto
		starts{0.0, -3.0, -10.0},
		ends{1.0, 2.0, 30.0};

	for (const auto start: starts) {
	for (const auto end: ends) {

		integral_type integral{0};
		const double step_size = (end - start) / 10;

		boost::numeric::odeint::integrate(
			// function to integrate
			[](const integral_type&, integral_type& dfdx, const double) {
				dfdx[0] = 3;
			},
			integral,
			start,
			end,
			step_size
		);

		const double
			i = integral[0],
			// exact solution
			e = 3 * (end - start);

		if (std::fabs(i - e) > 1e-9 * std::max(std::fabs(i), std::fabs(e))) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error when integrating from " << start
				<< " to " << end << ", exact: " << e << ", numerical: " << i
				<< std::endl;
			abort();
		}
	}}


	for (const auto start: starts) {
	for (const auto end: ends) {

		integral_type integral{0};
		const double step_size = (end - start) / 10;

		boost::numeric::odeint::integrate(
			[](const integral_type&, integral_type& dfdx, const double x) {
				dfdx[0] = std::pow(x, 3);
			},
			integral,
			start,
			end,
			step_size
		);

		const double
			i = integral[0],
			e = (std::pow(end, 4) - std::pow(start, 4)) / 4;

		if (std::fabs(i - e) > 1e-9 * std::max(std::fabs(i), std::fabs(e))) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error when integrating from " << start
				<< " to " << end << ", exact: " << e << ", numerical: " << i
				<< std::endl;
			abort();
		}
	}}


	for (const auto start: starts) {
	for (const auto end: ends) {

		integral_type integral{0};
		const double step_size = (end - start) / 10;

		boost::numeric::odeint::integrate(
			[](const integral_type&, integral_type& dfdx, const double x) {
				dfdx[0] = sin(x);
			},
			integral,
			start,
			end,
			step_size
		);

		const double
			i = integral[0],
			e = -cos(end) + cos(start);

		if (std::fabs(i - e) > 1e-5 * std::max(std::fabs(i), std::fabs(e))) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error when integrating from " << start
				<< " to " << end << ", exact: " << e << ", numerical: " << i
				<< std::endl;
			abort();
		}
	}}


	/* from https://en.wikipedia.org/wiki/Integral */
	integral_type integral{0};
	boost::numeric::odeint::integrate(
		[](const integral_type&, integral_type& dfdx, const double x) {
			dfdx[0]
				= (
					(322 + 3 * x * (98 + x * (37 + x))) / 100
					- 24 * x / (1 + x*x)
				) / 5;
		},
		integral,
		-2.0,
		2.0,
		1e-3
	);
	double
		i = integral[0],
		e = 94.0 / 25.0;
	if (std::fabs(i - e) > 1e-6 * std::max(std::fabs(i), std::fabs(e))) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Too large relative error when integrating from -2 to 2, exact: "
			<< e << ", numerical: " << i
			<< std::endl;
		abort();
	}


	// integrate f = x*y from -1,-1 to 2,2
	integral_type outer_integral{0};
	double current_x = 0;
	// https://stackoverflow.com/questions/22838052
	boost::numeric::odeint::integrate_const(
		boost::numeric::odeint::runge_kutta4<integral_type>(),
		[&](
			const integral_type&,
			integral_type& dfdx,
			const double
		) {
			integral_type inner_integral{0};

			boost::numeric::odeint::integrate(
				[&current_x](
					const integral_type&,
					integral_type& dfdy,
					const double y
				) {
					dfdy[0] = current_x * y;
				},
				inner_integral,
				-1.0,
				2.0,
				1e-3
			);

			dfdx[0] = inner_integral[0];
		},
		outer_integral,
		-1.0,
		2.0,
		1e-3,
		[&current_x](const integral_type&, const double x) {
			current_x = x;
		}
	);
	i = outer_integral[0],
	e = 2.25;
	if (std::fabs(i - e) > 1e-3 * std::max(std::fabs(i), std::fabs(e))) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Too large relative error when integrating from -1,-1 to 2,2; exact: " << e << ", numerical: " << i
			<< std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
