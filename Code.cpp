#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <iomanip>
using namespace std;

// function to make a k value
std::vector<double> find_k(double m1, double m2, double l1, double l2, double h, double v[4][1]) {
    
    // making the M matrix
    double M[2][2] = { { l1 * (m1 + m2),   m2 * l2 * cos(v[0][0] - v[1][0]) },   { l1 * cos(v[0][0] - v[1][0]),   l2 } } ;
    
    // finding the determinant
    double determinant = M[0][0] * M[1][1] - M[0][1] * M[1][0] ;
    
    // calculating what will be multiplied to 1/det
    double inv_M_temp[2][2] = { { M[1][1], - M[0][1] }, {- M[1][0], M[0][0]  } } ;
    
    
    // Making inv_M matrix
    double inv_M[2][2];
    
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            inv_M[i][j] = 1/determinant * inv_M_temp[i][j] ;
        }
    }
    
    
    // initialising the g vector
    double g[2][1] = { { -m2 * l2 * sin(v[0][0] - v[1][0]) * v[3][0]*v[3][0] - (m1 + m2) * 9.81 * sin(v[0][0]) }, { l1 * sin(v[0][0] - v[1][0]) * v[2][0]*v[2][0] - 9.81 * sin(v[1][0]) } } ;
    
    // Making the inv_M_times_g matrix
    double inv_M_times_g[2][1] = { {0.0}, {0.0} } ;
    
    // Filling in the inv_M_times_g matrix
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
            inv_M_times_g[j][0] = inv_M_times_g[j][0] + ( inv_M[j][i] * g[i][0] );
        }
    }
    
    double f[4][1] = { {v[2][0]}, {v[3][0]}, {inv_M_times_g[0][0]}, {inv_M_times_g[1][0]} } ;
    
    // Making k vector
    std::vector<double> k(4);
    
    // Filling k vector
    for (int i = 0; i < 4; ++i) {
        k[i] = h * f[i][0] ;
    }

    
    return k;
}


int main() {
    // opening the file
    ifstream vMyFile("parameters.txt");
    
    // making the parameters vector of size 8
    vector<double> parameters(8);
    
    // checking if the file is ok to read from
    if (vMyFile.good()) {
        
        int i = 0;
        
        while (true) {
            // reading the items in parameters
            vMyFile >> parameters[i];
            
            ++i;
            
            if (vMyFile.eof()) {
                break;
            }
            
        } 
        vMyFile.close();
    }
    else {
        cout << "Failed to open file" << endl;
    }
    
    
    // Doing Runge-Kutta
    
    // Making all the variables
    double m1, m2, l1, l2, v1, v2, h, t_f;
    
    // Making an array of pointers to the variables
    double *parameter_vars[8] = {&m1, &m2, &l1, &l2, &v1, &v2, &h, &t_f};

    
    // Making the variables equal to what is in the parameters vector
    for (int i = 0; i < parameters.size(); ++i) {        
        // Dereferencing the array of pointers so that the variables are given a value
        *parameter_vars[i] = parameters[i];
        
    }
    
    int size_of_y = t_f/h ;
    
    // initialising the y array
    double y[4][size_of_y];
    
    // the vector of initial conditions
    double v_0[4][1] = {{v1}, {v2}, {0}, {0}};
    
    double v[4][1];
    
    // Putting the initial values into the big y matrix
    for (int i = 0; i < 4; ++i) {
        y[i][0] = v_0[i][0];
    }
    
    
    // Opening output file
    ofstream vOutFile("output.txt", ios::out | ios::trunc);
    
    // Outputting the titles
    vOutFile.precision(4);
    vOutFile << setw(20) << "# time [s]"
    << setw(20) << "mass 1 (x, y) [m]"
    << setw(20) << "mass 2 (x, y) [m]" << endl;
    
    
    // Running the main loop of Runge-Kutta
    for (int main_i = 0; main_i <= t_f/h; ++main_i) {
        
        if (main_i < t_f/h) {
        
        // start of finding k_1
        
        // making v
        for (int i = 0; i < 4; ++i) {
            v[i][0] = y[i][main_i];
        }
        
        vector<double> k_1 = find_k(m1, m2, l1, l2, h, v) ;
        
        // end of finding k_1
        
        
        // start of finding k_2
        
        // making v
        for (int i = 0; i < 4; ++i) {
            v[i][0] = y[i][main_i] + 0.5 * k_1[i] ;
        }
        
        vector<double> k_2 = find_k(m1, m2, l1, l2, h, v) ;
        
        // end of finding k_2
        

        // start of finding k_3
        
        // making v
        for (int i = 0; i < 4; ++i) {
            v[i][0] = y[i][main_i] + 0.5 * k_2[i] ;
        }
        
        vector<double> k_3 = find_k(m1, m2, l1, l2, h, v) ;
        
        // end of finding k_3
        
        
        // start of finding k_4
        
        // making v
        for (int i = 0; i < 4; ++i) {
            v[i][0] = y[i][main_i] + k_3[i] ;
        }
        
        vector<double> k_4 = find_k(m1, m2, l1, l2, h, v) ;
        
        // end of finding k_4
        
        // Finding the next y value
        for (int i = 0; i < 4; ++i) {
            y[i][main_i + 1] = 0 ; // Setting y to 0 before changing it
            y[i][main_i + 1] = y[i][main_i] + 1.0/6.0 * k_1[i] + 1.0/3.0 * (k_2[i] + k_3[i]) + 1.0/6.0 * k_4[i];
            
        }
        
        }
        
        
        // Below happens for main_i up to and including t_f/h
        
        double time = h * (main_i) ;
        
        double x_1 = l1 * sin( y[0][main_i] ) ;
        double y_1 = - l1 * cos( y[0][main_i] ) ;
        
        double x_2 = l1 * sin( y[0][main_i] ) + l2 * sin( y[1][main_i] ) ;
        double y_2 = - l1 * cos( y[0][main_i] ) - l2 * cos( y[1][main_i] ) ;

        // Outputting data to output.txt
        vOutFile << setw(20) << time
        << setw(20) << "(" << x_1 << ", " << y_1 << ")"
        << setw(20) << "(" << x_2 << ", " << y_2 << ")" << endl;
       
    }
    
    // Closing the file
    vOutFile.close();
    
    return 0;
    
}


