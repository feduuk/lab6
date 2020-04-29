#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double func(double x1, double x2, double x3)
{
	return x1*x1 + 8 * x1 + 5 * x2*x2 + 7 * x3*x3 + 119.0*x3 + 531.75;
}

double f(double t, double x1, double x2, double x3, double grad1, double grad2, double grad3)
{
	return (x1 - t * grad1) * (x1 - t * grad1) + 8 * (x1 - t * grad1) + 5 * (x2 - t * grad2)*(x2 - t * grad2) + 7 * (x3 - t * grad3)*(x3 - t * grad3) + 119.0*(x3 - t * grad3) + 531.75;
}

double Method_kvadr_interpol(double a, double b, double h, double epsilon, double x1, double x2, double x3, double grad1, double grad2, double grad3)
{
	double t1 = (a + b) / 2;
	double t2;
	double t3;
	double num, denum;
	double t_star;
	double t_min;
	bool check = false;
	do
	{
		if (check)
		{
			if ((t_star > t1 && t_star < t3)  ||  (t_star > t3 && t_star < t1))
			{
				if (f(t_min, x1, x2, x3, grad1, grad2, grad3) < f(t_star, x1, x2, x3, grad1, grad2, grad3))
				{
					t1 = t_min;
				}
				else
				{
					t1 = t_star;
				}
			}
			else
			{
				t1 = t_star;
			}
		}

		t2 = t1 + h;
		if (f(t1, x1, x2, x3, grad1, grad2, grad3) > f(t2, x1, x2, x3, grad1, grad2, grad3))
		{
			t3 = t1 + 2 * h;
		}
		else
		{
			t3 = t1 - h;
		}

		if (f(t1, x1, x2, x3, grad1, grad2, grad3) < f(t2, x1, x2, x3, grad1, grad2, grad3))
		{
			if (f(t1, x1, x2, x3, grad1, grad2, grad3) < f(t3, x1, x2, x3, grad1, grad2, grad3))
			{
				t_min = t1;
			}
			else
			{
				t_min = t3;
			}
		}
		else
		{
			if (f(t2, x1, x2, x3, grad1, grad2, grad3) < f(t3, x1, x2, x3, grad1, grad2, grad3))
			{
				t_min = t2;
			}
			else
			{
				t_min = t3;
			}
		}

		num = (t2*t2 - t3 * t3) * f(t1, x1, x2, x3, grad1, grad2, grad3) + (t3*t3 - t1 * t1) * f(t2, x1, x2, x3, grad1, grad2, grad3) + (t1*t1 - t2 * t2) * f(t3, x1, x2, x3, grad1, grad2, grad3);
		denum = (t2 - t3) * f(t1, x1, x2, x3, grad1, grad2, grad3) + (t3 - t1) * f(t2, x1, x2, x3, grad1, grad2, grad3) + (t1 - t2) * f(t3, x1, x2, x3, grad1, grad2, grad3);

		if (denum == 0)
		{
			t1 = t_min;
			num = (t2*t2 - t3 * t3) * f(t1, x1, x2, x3, grad1, grad2, grad3) + (t3*t3 - t1 * t1) * f(t2, x1, x2, x3, grad1, grad2, grad3) + (t1*t1 - t2 * t2) * f(t3, x1, x2, x3, grad1, grad2, grad3);
			denum = (t2 - t3) * f(t1, x1, x2, x3, grad1, grad2, grad3) + (t3 - t1) * f(t2, x1, x2, x3, grad1, grad2, grad3) + (t1 - t2) * f(t3, x1, x2, x3, grad1, grad2, grad3);
		}

		t_star = 0.5 * (num / denum);
		check = true;
	} while ((abs(f(t_min, x1, x2, x3, grad1, grad2, grad3) - f(t_star, x1, x2, x3, grad1, grad2, grad3)) > epsilon) || (abs(t_min - t_star) > epsilon));

	return t_star;
}

void Skolz_okno(double &a, double &b, double x0, double h, double x1, double x2, double x3, double grad1, double grad2, double grad3)
{
	if ((f(x0, x1, x2, x3, grad1, grad2, grad3) < f((x0 - h), x1, x2, x3, grad1, grad2, grad3)) && (f(x0, x1, x2, x3, grad1, grad2, grad3) < f((x0 + h), x1, x2, x3, grad1, grad2, grad3)))
	{
		a = x0 - h;
		b = x0 + h;
	}
	else if(  f(x0 - h, x1, x2, x3, grad1, grad2, grad3) > f((x0 + h), x1, x2, x3, grad1, grad2, grad3)  )
	{
		x0 = x0 + h;
		Skolz_okno(a, b, x0, h, x1, x2, x3, grad1, grad2, grad3);
	}
	else
	{
		x0 = x0 - h;
		Skolz_okno(a, b, x0, h, x1, x2, x3, grad1, grad2, grad3);
	}
}

void Method_naiskor_spuska(double x01, double x02, double x03, double &x1_finish, double &x2_finish, double &x3_finish)
{
	double x1 = x01;
	double x2 = x02;
	double x3 = x03;
	double x1_new, x2_new, x3_new;
	double t;
	double grad1, grad2, grad3;
	double a;
	double b;
	double epsilon = 0.01;
	double h = 0.01;
	double x0 = 0.5;
	bool check = false;
	
	do{
		if (check)
		{
			x1 = x1_new;
			x2 = x2_new;
			x3 = x3_new;
		}
		grad1 = 2 * x1 + 8;
		grad2 = 10 * x2;
		grad3 = 14 * x3 + 119.0;
		Skolz_okno(a, b, x0, h, x1, x2, x3, grad1, grad2, grad3);
		t = Method_kvadr_interpol(a, b, h, epsilon, x1, x2, x3, grad1, grad2, grad3);
		x1_new = x1 - t * grad1;
		x2_new = x2 - t * grad2;
		x3_new = x3 - t * grad3;
		check = true;
	} while ((abs(func(x1_new, x2_new, x3_new) - func(x1, x2, x3)) > epsilon) || (abs(x1_new - x1) > epsilon) || (abs(x2_new - x2) > epsilon) || (abs(x3_new - x3) > epsilon));

	x1_finish = x1_new;
	x2_finish = x2_new;
	x3_finish = x3_new;
}
int main()
{
	ofstream fout1("C:\\Users\\Федор\\Desktop\\file1.txt");
	double x1, x2, x3;
	Method_naiskor_spuska(0, 0, 0, x1, x2, x3);
	cout <<"x1 = " << x1 << "    x2 = " << x2 << "    x3 = " << x3 << endl;
	cout << "f(x1,x2,x3) = " << func(x1, x2, x3) << endl;

	//double x11, x22, x33;
	//double foo = 100;
	//x33 = x3;
	//x11*x11 + 8 * x11 + 5 * x22*x22 + 7 * x33*x33 + 119.0*x33 + 531.75;

	//while (foo < 1000)
	//{
	//	for (int i = 0; i < 6; i++)
	//	{
	//		x11 = i;
	//		x22 = pow((-x11 * x11 - 8 * x11 - 7 * x33*x33 - 119.0*x33 - 531.75 + foo) / 5, 0.5);
	//		fout1 << x22 << "   " << x11 << endl;
	//	}
	//	foo += 100;
	//	break;
	//}
	fout1.close();
	system("pause");
	return 0;
}