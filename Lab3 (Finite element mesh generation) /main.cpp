#include <iostream>
#include <vector>
#include <fstream>
#include <string>

struct Point
{
	double x = 0;
	double y = 0;
	double z = 0;
	int number = 0;
};

struct Element
{
	std::vector<Point> points;
	int number = 0;
};

std::vector<double> gen_grid(const int& m, const double& t0, const double& tn)
{
	double h = (tn - t0) / (m-1);
	std::vector<double> z(m);
	for (int i = 0; i < m; i++)
	{
		z[i] = t0 + i * h;
	}
	return z;
}


std::vector<Point> linear_type1(const std::vector<Point>& coordinates, int& m)
{
	m++;
	std::vector<Point> new_coordinates;
	std::vector<double> x = gen_grid(m, coordinates[0].x, coordinates[1].x);
	std::vector<double> y = gen_grid(m, coordinates[0].y, coordinates[1].y);
	std::vector<double> z = gen_grid(m, coordinates[0].z, coordinates[1].z);
	for (int i = 0; i < m; i++)
	{
		new_coordinates.emplace_back(Point({x[i],y[i],z[i],i+1}));
	}
	return new_coordinates;
}

std::vector<Point> linear_type2(const std::vector<Point>& coordinates, int& m)
{
	int size = 1;
	m++;
	std::vector<Point> new_coordinates;
	std::vector<double> x = gen_grid(m, coordinates[0].x, coordinates[1].x);
	std::vector<double> y = gen_grid(m, coordinates[0].y, coordinates[1].y);
	std::vector<double> z = gen_grid(m, coordinates[0].z, coordinates[1].z);
	for (int i = 0; i < m-1; i++)
	{
		new_coordinates.emplace_back(Point({x[i], y[i], z[i], i + size}));
		size++;
		new_coordinates.emplace_back(Point({(x[i] + x[i + 1]) / 2, (y[i] + y[i + 1]) / 2, (z[i] + z[i + 1]) / 2, i + size}));
	}
	new_coordinates.emplace_back(Point({x[m - 1], y[m - 1], z[m - 1], 2 * m - 1}));
	return new_coordinates;
}

std::vector<Point> quadrilateral(const std::vector<Point>& coordinates, int& m, int&n)
{
	// n - по ветикали, m - по горизонтали
	m++;
	n++;
	std::vector<Point> new_coordinates;
	std::vector<double> x_left = gen_grid(n, coordinates[0].x, coordinates[3].x);
	std::vector<double> y_left = gen_grid(n, coordinates[0].y, coordinates[3].y);
	std::vector<double> z_left = gen_grid(n, coordinates[0].z, coordinates[3].z);

	std::vector<double> x_right = gen_grid(n, coordinates[1].x, coordinates[2].x);
	std::vector<double> y_right = gen_grid(n, coordinates[1].y, coordinates[2].y);
	std::vector<double> z_right = gen_grid(n, coordinates[1].z, coordinates[2].z);

	std::vector<double> x_tmp;
	std::vector<double> y_tmp;
	std::vector<double> z_tmp;
	for (int i = 0; i < n; i++)
	{
		x_tmp = gen_grid(m, x_left[i], x_right[i]);
		y_tmp = gen_grid(m, y_left[i], y_right[i]);
		z_tmp = gen_grid(m, z_left[i], z_right[i]);
		for (int j = 0; j < m; j++)
		{
			new_coordinates.emplace_back(Point({x_tmp[j], y_tmp[j], z_tmp[j], i * m + j + 1}));
		}

	}
	return new_coordinates;
}

std::vector<Point> quadrilateral_type4(const std::vector<Point>& coordinates, int& m, int& n)
{
	// n - по ветикали, m - по горизонтали
	m++;
	n++;
	std::vector<Point> new_coordinates;
	std::vector<double> x_left = gen_grid(n, coordinates[0].x, coordinates[3].x);
	std::vector<double> y_left = gen_grid(n, coordinates[0].y, coordinates[3].y);
	std::vector<double> z_left = gen_grid(n, coordinates[0].z, coordinates[3].z);

	std::vector<double> x_right = gen_grid(n, coordinates[1].x, coordinates[2].x);
	std::vector<double> y_right = gen_grid(n, coordinates[1].y, coordinates[2].y);
	std::vector<double> z_right = gen_grid(n, coordinates[1].z, coordinates[2].z);

	std::vector<double> x_tmp;
	std::vector<double> y_tmp;
	std::vector<double> z_tmp;
	for (int i = 0; i < n; i++)
	{
		x_tmp = gen_grid(m, x_left[i], x_right[i]);
		y_tmp = gen_grid(m, y_left[i], y_right[i]);
		z_tmp = gen_grid(m, z_left[i], z_right[i]);
		for (int j = 0; j < m; j++)
		{
			new_coordinates.emplace_back(Point({x_tmp[j], y_tmp[j], z_tmp[j], i * m + j + 1 }));
		}

	}
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	int number = new_coordinates.size();
	for (int i = 0; i < (n-1) * (m-1); i++)
	{
		sum1 = new_coordinates[i + (i / (m - 1))].x + new_coordinates[i + 1 + (i / (m - 1))].x + new_coordinates[i + 1 + m + (i / (m - 1))].x + new_coordinates[i + m + (i / (m - 1))].x;
		sum2 = new_coordinates[i + (i / (m - 1))].y + new_coordinates[i + 1 + (i / (m - 1))].y + new_coordinates[i + 1 + m + (i / (m - 1))].y + new_coordinates[i + m + (i / (m - 1))].y;
		sum3 = new_coordinates[i + (i / (m - 1))].z + new_coordinates[i + 1 + (i / (m - 1))].z + new_coordinates[i + 1 + m + (i / (m - 1))].z + new_coordinates[i + m + (i / (m - 1))].z;
		new_coordinates.emplace_back(Point({sum1 / 4, sum2 / 4, sum3 / 4, number + i + 1}));
	}
	return new_coordinates;
}


std::vector<Element> get_elements(const int& num_points, const int& type, const int& NE1, const int& NE2, const std::vector<Point>& coordinates)
{
	std::vector<Element> elements;
	int size = 0;
	if (num_points == 4){
		switch (type) {
			case 1:
				elements.resize((NE1 - 1) * (NE2 - 1));
				for (int i = 0; i < (NE1 - 1) * (NE2 - 1); i++) {
					elements[i].number = i + 1;
					elements[i].points.emplace_back(coordinates[i + (i / (NE1 - 1))]);
					elements[i].points.emplace_back(coordinates[i + 1 + (i / (NE1 - 1))]);
					elements[i].points.emplace_back(coordinates[i + 1 + NE1 + (i / (NE1 - 1))]);
					elements[i].points.emplace_back(coordinates[i + NE1 + (i / (NE1 - 1))]);
				}
				break;
			case 2:
				elements.resize(2*(NE1 - 1) * (NE2 - 1));
				for (int i = 0; i < (NE1 - 1) * (NE2 - 1); i++) {
					elements[size].number = i*2 + 1;
					elements[size].points.emplace_back(coordinates[i + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + 1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + 1 + NE1 + (i / (NE1 - 1))]);
					size++;
					elements[size].number = i * 2 + 2;
					elements[size].points.emplace_back(coordinates[i + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + 1 + NE1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + NE1 + (i / (NE1 - 1))]);
					size++;
				}
				break;
			case 3:
				elements.resize(2 * (NE1 - 1) * (NE2 - 1));
				for (int i = 0; i < (NE1 - 1) * (NE2 - 1); i++) {
					elements[size].number = i * 2 + 1;
					elements[size].points.emplace_back(coordinates[i + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + 1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + NE1 + (i / (NE1 - 1))]);
					size++;
					elements[size].number = i * 2 + 2;
					elements[size].points.emplace_back(coordinates[i + 1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + 1 + NE1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + NE1 + (i / (NE1 - 1))]);
					size++;
				}
				break;
			case 4:
				elements.resize(4 * (NE1 - 1) * (NE2 - 1));
				for (int i = 0; i < (NE1 - 1) * (NE2 - 1); i++) {
					elements[size].number = i * 4 + 1;
					elements[size].points.emplace_back(coordinates[i + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + 1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[NE1*NE2 + i]);
					size++;
					elements[size].number = i * 4 + 2;
					elements[size].points.emplace_back(coordinates[i + 1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + 1 + NE1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[NE1 * NE2 + i]);
					size++;
					elements[size].number = i * 4 + 3;
					elements[size].points.emplace_back(coordinates[i + 1 + NE1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + NE1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[NE1 * NE2 + i]);
					size++;
					elements[size].number = i * 4 + 4;
					elements[size].points.emplace_back(coordinates[i + NE1 + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[i + (i / (NE1 - 1))]);
					elements[size].points.emplace_back(coordinates[NE1 * NE2 + i]);
					size++;
				}
				break;
		}
	}
	else
	{
		elements.resize(NE1 - 1);
		if (type == 1){
			for (int i = 0; i < NE1 - 1; i++){
				for (int j = 0; j < 2; j++){
					elements[i].points.emplace_back(coordinates[j+i]);
					elements[i].number = i + 1;
				}
			}
		}
		else{
			for (int i = 0; i < NE1 - 1; i++){
				for (int j = 0; j < 3; j++){
					elements[i].points.emplace_back(coordinates[j + i*2]);
					elements[i].number = i + 1;
				}
			}
		}
	}
	return elements;
}

void generate_grid(int num_points, int type, int& NE1, int& NE2, std::vector<Point>& initial_coordinates, std::vector<Point>& coordinates, std::vector<Element>& FE)
{
 	if (num_points == 4) {
		if (type == 4)
		{
			coordinates = quadrilateral_type4(initial_coordinates, NE1, NE2);
		}
		else
		{
			coordinates = quadrilateral(initial_coordinates, NE1, NE2);
		}
	}
	else
	{
		if (type == 1)
		{
			coordinates = linear_type1(initial_coordinates, NE1);
		}
		else
		{
			coordinates = linear_type2(initial_coordinates, NE1);
		}
	}
	FE = get_elements(num_points, type, NE1, NE2, coordinates);
}

std::vector<int> get_contour(const int& NE1, const int& NE2,const std::vector<Point>& coordinates)
{
	std::vector<int> contour_numbers;
	for (int i = 0; i < NE1; i++){contour_numbers.emplace_back(coordinates[i].number);}
	for (int i = 1; i < NE2; i++) {contour_numbers.emplace_back(coordinates[(i + 1) * (NE1) - 1].number); }
	for (int i = 1; i < NE1; i++) {contour_numbers.emplace_back(coordinates[(NE2) * (NE1) - 1 - i].number); }
	for (int i = 1; i < NE2-1; i++) {contour_numbers.emplace_back(coordinates[NE1 * i].number); }
	return contour_numbers;
}

void export_to_file(const std::string& file_output, const int& num_points ,const int& type, const int& NE1, const int& NE2, std::vector<Point>& coordinates, std::vector<Element>& FE)
{
	std::ofstream out(file_output);
	if (out.is_open()) {
		if (num_points == 4){
			out << FE.size() << ' ';
			switch (type){
				case 1:
					out << (NE1) * (NE2) << ' ';
					break;
				case 2:
					out << (NE1) * (NE2) << ' ';
					break;
				case 3:
					out << (NE1) * (NE2) << ' ';
					break;
				case 4:
					out << (NE1) * (NE2) + (NE1 - 1)*(NE2-1)<< ' ';
					break;
			}
			out << 1 << '\n';
		}
		else
		{
			out << FE.size() << ' ';
			if (type == 1)
			{
				out << NE1 <<' ';
			}
			else
			{
				out << 2*NE1 - 1 << ' ';
			}
			out << 2 << '\n';
		}
		//выводим номера элементов
		for (int i = 0; i < FE.size(); i++)
		{
			out << FE[i].number << ' ' ;
			out << FE[i].points.size() << ' ' ;
			for (int j = 0; j < FE[i].points.size() - 1; j++){out << FE[i].points[j].number << ' ';}
			out << FE[i].points[FE[i].points.size() - 1].number << '\n';
		}

		for (int i = 0; i < coordinates.size(); i++)
		{
			out << coordinates[i].number << ' ';
			out << coordinates[i].x << ' ';
			out << coordinates[i].y << ' ';
			out << coordinates[i].z << '\n';
		}
		if (num_points == 4) {
			auto CountourNum = get_contour(NE1, NE2, coordinates);
			out << CountourNum.size() << '\n';
			for (int i = 0; i < CountourNum.size(); i++)
			{
				out << CountourNum[i] << '\n';
			}
		}
		else {
			out << 1 << ' ' << 1 << '\n';
			out << coordinates[0].number << '\n';
			out << coordinates[coordinates.size() - 1].number << '\n';
		}

	}
	out.close();
}
int main()
{
	int num_points = 0;
	int type = 0;
	int NE1 = 0;
	int NE2 = 0;

	std::string input_file("test2.txt");
	std::ifstream file(input_file);
	file >> num_points;
	std::vector<Point> coordinates(num_points);
	for (int i = 0; i < num_points; i++)
	{
		file >> coordinates[i].x >> coordinates[i].y >> coordinates[i].z;
	}
	if (num_points == 4)
	{
		file >> NE1 >> NE2;
	}
	else
	{
		file >> NE1;
	}
	file >> type;
	file.close();

	std::vector<Point> fcoordinates;
	std::vector<Element> FE;
	generate_grid(num_points, type, NE1, NE2, coordinates, fcoordinates, FE);
	export_to_file("output2.txt", num_points, type, NE1, NE2, fcoordinates, FE);
}