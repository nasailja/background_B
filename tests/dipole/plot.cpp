/*
Plots a few tests for dipole code.

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
#include "cstdlib"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "fstream"
#include "iostream"
#include "utility"
#include "vector"

#include "dipole.hpp"
#include "volume_range.hpp"

using namespace std;
using namespace background_B;
using namespace pamhd;

int main()
{
	// radius of earth
	constexpr double Re = 6.3712e6;
	constexpr std::array<double, 2>
		sample_start{-10 * Re, -10 * Re},
		sample_end{10 * Re, 10 * Re};
	constexpr size_t samples = 9;

	std::vector<
		std::pair<
			Eigen::Vector3d,
			Eigen::Vector3d
		>
	> dipole_moments_positions{{
		get_earth_dipole_moment<Eigen::Vector3d>(),
		{0, 0, 0}
	}};

	// plot output with gnuplot

	// plot 1
	ofstream plot1("plot_y=0_dip=0,0,0.gnuplot", ios_base::out);
	plot1 <<
		"set term svg enhanced\n"
		"set output 'plot_y=0_dip=0,0,0.svg'\n"
		"set size square\n"
		"plot '-' using 1:2:3:4 with vectors title ''\n";
	for (const auto& coord: volume_range(sample_start, sample_end, samples)) {
		const Eigen::Vector3d
			field = get_dipole_field_0d(
				dipole_moments_positions,
				{coord[0], 0, coord[1]}
			),
			nomalized_field = field / field.norm();

		plot1 << coord[0] / Re << " " << coord[1] / Re << " "
			<< nomalized_field[0] << " " << nomalized_field[2] << "\n";
	}
	plot1 << "end" << endl;
	plot1.close();
	system("gnuplot plot_y=0_dip=0,0,0.gnuplot");

	// plot 2
	dipole_moments_positions[0].second << 10 * Re, 1 * Re, 10 * Re;
	ofstream plot2("plot_y=0_dip=10,1,10.gnuplot", ios_base::out);
	plot2 <<
		"set term svg enhanced\n"
		"set output 'plot_y=0_dip=10,1,10.svg'\n"
		"set size square\n"
		"plot '-' using 1:2:3:4 with vectors title ''\n";
	for (const auto& coord: volume_range(sample_start, sample_end, samples)) {
		const Eigen::Vector3d
			field = get_dipole_field_0d(
				dipole_moments_positions,
				{coord[0], 0, coord[1]}
			),
			nomalized_field = field / field.norm();

		plot2 << coord[0] / Re << " " << coord[1] / Re << " "
			<< nomalized_field[0] << " " << nomalized_field[2] << "\n";
	}
	plot2 << "end" << endl;
	plot2.close();
	system("gnuplot plot_y=0_dip=10,1,10.gnuplot");

	// plot 3
	dipole_moments_positions[0].second << 0, 0, 0;
	const Eigen::Matrix3d rotator(
		Eigen::AngleAxisd(1, Eigen::Vector3d::UnitY())
		* Eigen::AngleAxisd(1, Eigen::Vector3d::UnitZ())
	);
	dipole_moments_positions[0].first
		= (rotator * dipole_moments_positions[0].first).eval();
	ofstream plot3("plot_y=0_rot.gnuplot", ios_base::out);
	plot3 <<
		"set term svg enhanced\n"
		"set output 'plot_y=0_rot.svg'\n"
		"set size square\n"
		"plot '-' using 1:2:3:4 with vectors title ''\n";
	for (const auto& coord: volume_range(sample_start, sample_end, samples)) {
		const Eigen::Vector3d
			field = get_dipole_field_0d(
				dipole_moments_positions,
				{coord[0], 0, coord[1]}
			),
			nomalized_field = field / field.norm();

		plot3 << coord[0] / Re << " " << coord[1] / Re << " "
			<< nomalized_field[0] << " " << nomalized_field[2] << "\n";
	}
	plot3 << "end" << endl;
	plot3.close();
	system("gnuplot plot_y=0_rot.gnuplot");

	// plot 4
	dipole_moments_positions[0].second << -7 * Re, 0, 0;
	dipole_moments_positions.push_back({
		get_earth_dipole_moment<Eigen::Vector3d>(),
		{7 * Re, 0, 0}
	});
	ofstream plot4("plot_y=0_rot_2.gnuplot", ios_base::out);
	plot4 <<
		"set term svg enhanced\n"
		"set output 'plot_y=0_rot_2.svg'\n"
		"set size square\n"
		"set xrange [-11 : 11]\n"
		"set yrange [-11 : 11]\n"
		"plot '-' using 1:2:3:4 with vectors title ''\n";
	for (const auto& coord: volume_range(sample_start, sample_end, 2 * samples)) {
		const Eigen::Vector3d
			field = get_dipole_field_0d(
				dipole_moments_positions,
				{coord[0], 0, coord[1]}
			),
			nomalized_field = field / field.norm();

		plot4 << coord[0] / Re << " " << coord[1] / Re << " "
			<< nomalized_field[0] << " " << nomalized_field[2] << "\n";
	}
	plot4 << "end" << endl;
	plot4.close();
	system("gnuplot plot_y=0_rot_2.gnuplot");

	return EXIT_SUCCESS;
}
