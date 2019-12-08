#include <stdio.h>
#include <math.h>

/*
 ( x')          (x)
 ( y')  =   M * (y)
i ( 1 )          (1)

instead of (x',y',1) = (x,y,1) * M
*/

// here's how to print a 3x3 array
int M2d_print_mat (double a[3][3]) {
  for (int r = 0; r < 3; r++) {
      for (int c = 0; c < 3; c++) {
           printf(" %12.4lf ", a[r][c]);
      }
      printf("\n");
  }

  return 1;
}

// copy matrix
int M2d_copy_mat (double a[3][3], double b[3][3]) {
    // a = b
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = b[i][j];
        }
    }

    return 1;
}

// all spots are 0, unless row = column
int M2d_make_identity (double a[3][3]) {
    // a = I
    int r, c;
    for (r = 0; r < 3; r++) {
        for (c = 0; c < 3; c++) {
            if (r == c) {
                a[r][c] = 1.0;
            }
            else {
                a[r][c] = 0.0;
            }
        }
    }

    return 1;
}

// do a translation
int M2d_make_translation (double a[3][3], double dx, double dy) {
    M2d_make_identity(a) ;

    a[0][2] = dx;
    a[1][2] = dy;

    return 1;
}

// multiply two 3x3 matrices
int M2d_mat_mult (double res[3][3], double a[3][3], double b[3][3]) {
    // res = a * b
    // this is SAFE, i.e. the user can make a call such as
    // M2d_mat_mult(p, p, q) or M2d_mat_mult(p, q, p) or  M2d_mat_mult(p, p, p)

    double first[3][3];
    double second[3][3];

    M2d_copy_mat(first, a);
    M2d_copy_mat(second,b);


    // math
    for (int first_row = 0; first_row < 3; first_row++) {
        for (int second_col = 0; second_col < 3; second_col++) {
            double result_point = 0;

            for (int point = 0; point < 3; point++) {
                result_point += first[first_row][point] * second[point][second_col];
            }

            res[first_row][second_col] = result_point;
        }
    }

    return 1;
}

// make a scale
int M2d_make_scaling (double a[3][3], double sx, double sy) {
    M2d_make_identity(a);
    a[0][0] = sx;
    a[1][1] = sy;
    // a[2][2] = 1; will happen in identity

    return 1;
}

// make rotation in radians
int M2d_make_rotation_radians (double a[3][3],  double radians) {
    // rotation matrix
    // [cos(radians)   -sin(radians)]
    // [sin(radians)    cos(radians)]
    M2d_make_identity(a);
    a[0][0] = cos(radians);
    a[0][1] = -sin(radians);
    a[1][0] = sin(radians);
    a[1][1] = cos(radians);

    return 1;
}

// make rotation in degrees
int M2d_make_rotation_degrees (double a[3][3],  double degrees) {
    M2d_make_rotation_radians(a, degrees * (M_PI / 180));

    return 1;
}

// make rotation, assuming that you already know sine and cosine of the values
int M2d_make_rotation_cs (double a[3][3], double cs, double sn) {
    M2d_make_identity(a);
    a[0][0] = cs;
    a[0][1] = -sn;
    a[1][0] = sn;
    a[1][1] = cs;

    return 1;
}

int M2d_mat_mult_pt (double P[2], double m[3][3], double Q[2]) {
    // P = m*Q
    // SAFE, user may make a call like M2d_mat_mult_pt (W, m,W) ;

    // assume that the extra spot is a 1
    // [1, 2, 3]   [10]
    // [4, 5, 6] * [11]
    // [7, 8, 9]   [ 1]

    double point[2];
    point[0] = Q[0];
    point[1] = Q[1];


    P[0] = m[0][0] * point[0] + m[0][1] * point[1] + m[0][2] * 1;
    P[1] = m[1][0] * point[0] + m[1][1] * point[1] + m[1][2] * 1;

    return 1;
}


int M2d_mat_mult_points (double *X, double *Y, double m[3][3], double *x, double *y, int numpoints) {
    // |X0 X1 X2 ...|       |x0 x1 x2 ...|
    // |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
    // | 1  1  1 ...|       | 1  1  1 ...|
    // SAFE, user may make a call like M2d_mat_mult_points (x,y, m, x,y, n) ;
  double temp ;

    for (int i = 0; i < numpoints; i++) {
      temp = m[0][0]*x[i] + m[0][1]*y[i] + m[0][2] ;
      Y[i] = m[1][0]*x[i] + m[1][1]*y[i] + m[1][2] ;
      X[i] = temp ;
    }

    return 1;
}

/*int main() {
    double a[3][3], b[3][3], result[3][3];
    a[0][0] = 1;
    a[0][1] = 2;
    a[0][2] = 3;
    a[1][0] = 4;
    a[1][1] = 5;
    a[1][2] = 6;
    a[2][0] = 1;
    a[2][1] = 1;
    a[2][2] = 1;

    b[0][0] = 7;
    b[0][1] = 8;
    b[0][2] = 2;
    b[1][0] = 9;
    b[1][1] = 10;
    b[1][2] = 2;
    b[2][0] = 11;
    b[2][1] = 12;
    b[2][2] = 2;

    M2d_mat_mult(result, a, b);

    M2d_print_mat(result);

    printf("\n");

    M2d_make_rotation_cs(a, 1, 0);

    M2d_print_mat(a);
}*/
