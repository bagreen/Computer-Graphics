// don't draw backfaces
// find vector perp to polygon
// find vector from polygon to viewer
// compute perp vector * viewer vector

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <FPT.h>
#include "M3d_matrix_toolsS.c"

int window_size = 800;
int num_points[10000], num_polys[10000];
int psize[50][1000];
int con[100][1000][1000];
double x[10][50000], y[10][50000], z[10][50000];

_Bool reverse_back  = False;

double dot_product(double vector1[3], double vector2[3]) {
    return ((vector1[0] * vector2[0]) + (vector1[1] * vector2[1]) + (vector1[2] * vector2[2]));
}

int cross_product(double result[3], double vector1[3], double vector2[3]) {
    result[0] =   (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
    result[1] = -((vector1[0] * vector2[2]) - (vector1[2] * vector2[0]));
    result[2] =   (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);
}

int make_vector(double vector[3], double x1, double y1, double z1, double x2, double y2, double z2) {
    vector[0] = x2 - x1;
    vector[1] = y2 - y1;
    vector[2] = z2 - z1;
}

int make_perp_vector(double perp_vector[3],
		     double x1, double y1, double z1,
		     double x2, double y2, double z2,
		     double x3, double y3, double z3) {
    double vector1[3], vector2[3];
    make_vector(vector1,  x1, y1, z1, x2, y2, z2);
    make_vector(vector2, x1, y1, z1, x3, y3, z3);
    
    cross_product(perp_vector, vector1, vector2);
}

int find_poly_center(double center[3], double *x, double *y, double *z, int length) {
    double x_total = 0, y_total = 0, z_total = 0;

    for (int i = 0; i < length; i++) {
        x_total += x[i];
        y_total += y[i];
        z_total += z[i];
    }

    center[0] = x_total / length;
    center[1] = y_total / length;
    center[2] = z_total / length;
}

int find_shape_center(int file_num) {
    // would just looping through all of the points in the polygon work as well/better? Or is this better? This can count points multiple times as they're used
    // find center through center of mass
    double x_total = 0, y_total = 0, z_total = 0;
    int count = 0;

    // get the x, y values of each point that is in each polygon
    for (int i = 0; i < num_polys[file_num]; i++) {
        for (int j = 0; j < psize[file_num][i]; j++) {
            x_total += x[file_num][con[file_num][i][j]];
            y_total += y[file_num][con[file_num][i][j]];
            z_total += z[file_num][con[file_num][i][j]];
            count   += 1;
        }
    }

    x[file_num][num_points[file_num]] = x_total / count;
    y[file_num][num_points[file_num]] = y_total / count;
    z[file_num][num_points[file_num]] = z_total / count;
}

void draw_polygon(int file_num) {
    // draw each polygon in the file
    for (int poly_num = 0; poly_num < num_polys[file_num]; poly_num++) {
        double x_points[1000], y_points[1000], z_points[1000];

        // get the x, y values of each point that is in each polygon
        for (int point = 0; point < psize[file_num][poly_num]; point++) {
            x_points[point] = x[file_num][con[file_num][poly_num][point]];
            y_points[point] = y[file_num][con[file_num][poly_num][point]];
            z_points[point] = z[file_num][con[file_num][poly_num][point]];
        }

        // find center of polygon
        double center[3];
        find_poly_center(center, x_points, y_points, z_points, psize[file_num][poly_num]);

        // make viewer vector
        double view_vector[3];
        if (reverse_back == True) {
            make_vector(view_vector, 0, 0, 0, center[0], center[1], center[2]);
        }
        else {
            make_vector(view_vector, center[0], center[1], center[2], 0, 0, 0);
        }

        // make perpindicular vector
        double perp_vector[3], first_pt[3], second_pt[3], third_pt[3];
        make_perp_vector(perp_vector,
			 x_points[0], y_points[0], z_points[0],
			 x_points[1], y_points[1], z_points[1],
			 x_points[2], y_points[2], z_points[2]); 

        // get dot product, if result is equal to/greater than 0, draw!
        double dot = dot_product(perp_vector, view_vector);
        if (dot > 0) {
            // do math to calculate 3d aspect of this shape
            double x_bar, y_bar, x_barbar, y_barbar, h;
            h = tan(45*M_PI/180);

            for (int point = 0; point < psize[file_num][poly_num]; point++) {
                x_bar = x_points[point] / z_points[point];
                y_bar = y_points[point] / z_points[point];

                x_barbar = (window_size / 2) * (x_bar / h) + (window_size / 2);
                y_barbar = (window_size / 2) * (y_bar / h) + (window_size / 2);

                x_points[point] = x_barbar;
                y_points[point] = y_barbar;
            }

            // print polygon
            G_polygon(x_points, y_points, psize[file_num][poly_num]);
        }
    }
}

int main(int argc, char **argv) {
	// get arguments on command line and check for files existence
    FILE *fp;
	int num_objects = argc - 1;

    // setup polygons, points, and numbers for each on the command line
   	for (int k = 0; k < num_objects; k++) {
		fp = fopen(argv[k + 1], "r");

        // if file doesn't exist, print error and quit
		if (fp == NULL) {
			printf("Can't open file, %s\n", argv[2]);
			exit(0);
		}

        // else, scan in values
        else {
		    // first line is the number of points
		    fscanf(fp, "%d", &num_points[k]);

            // next num_points lines are the x and y values of each point
		    for (int i = 0; i < num_points[k]; i++) {
                fscanf(fp, "%lf %lf %lf", &x[k][i], &y[k][i], &z[k][i]);
            }

            // next line is the number of polygons
		    fscanf(fp, "%d", &num_polys[k]);

            // next num_poly lines are the number of lines in this polygon
		    for (int i = 0; i < num_polys[k]; i++) {
			    fscanf(fp, "%d", &psize[k][i]);

                // followed by the points that are connected in this polygon
			    for (int j = 0; j < psize[k][i]; j++) {
                    fscanf(fp, "%d", &con[k][i][j]);
                }
		    }
        }
	}

    // make movement matrices
    // translate 1 in each direction matrices
    double x_trans_forw[4][4], y_trans_forw[4][4], z_trans_forw[4][4];
    double x_trans_back[4][4], y_trans_back[4][4], z_trans_back[4][4];
    M3d_make_translation(x_trans_forw,  0.1,  0,  0);
    M3d_make_translation(y_trans_forw,  0,  0.1,  0);
    M3d_make_translation(z_trans_forw,  0,  0,  0.1);
    M3d_make_translation(x_trans_back, -0.1,  0,  0);
    M3d_make_translation(y_trans_back,  0, -0.1,  0);
    M3d_make_translation(z_trans_back,  0,  0, -0.1);

    // rotate two degrees in each direction matrices
    double sine = sin((2 * M_PI) / 180);
    double cose = cos((2 * M_PI) / 180);
    double x_rot_forw[4][4], y_rot_forw[4][4], z_rot_forw[4][4];
    double x_rot_back[4][4], y_rot_back[4][4], z_rot_back[4][4];
    M3d_make_rotation_x(x_rot_forw,  sine, cose);
    M3d_make_rotation_y(y_rot_forw,  sine, cose);
    M3d_make_rotation_z(z_rot_forw,  sine, cose);
    M3d_make_rotation_x(x_rot_back, -sine, cose);
    M3d_make_rotation_y(y_rot_back, -sine, cose);
    M3d_make_rotation_z(z_rot_back, -sine, cose);

	// setup graphics, draw black screen
    G_init_graphics(window_size, window_size);
    G_rgb(0, 0, 0);
    G_fill_rectangle(0, 0, window_size, window_size);

    // translate shapes to origin
    for (int file_num = 0; file_num < num_objects; file_num++) {
        // storing center in x, y, z, so that operations on the figure will also move the center
        find_shape_center(file_num);
        double translate_origin[4][4];

        M3d_make_translation(translate_origin, -x[file_num][num_points[file_num]], -y[file_num][num_points[file_num]], -z[file_num][num_points[file_num]]);
        M3d_mat_mult_points(x[file_num], y[file_num], z[file_num], translate_origin, x[file_num], y[file_num], z[file_num], num_points[file_num] + 1);
    }

    int file_num = 0;
    int last_num = 0;

    _Bool reverse   = False;
    _Bool rotate    = False;
    _Bool translate = False;

    draw_polygon(file_num);

    // loop for getting keys and drawing polygons
    while (1) {
        double transform[4][4];
        M3d_make_identity(transform);

        int key = G_wait_key();
        int key_num = key - 49;

        printf("Key_num is %d\n", key_num);

        // display new object
        if (key_num < num_objects) {
            file_num = key_num;
        }

        // reverse perp vector for image
        else if (key == 'b') {
            if (reverse_back == True)    reverse_back = False;
            else                        reverse_back = True;
        }

        // reverse the transform/rotation
        else if (key == 'c') {
            if (reverse == True)        reverse = False;
            else                        reverse = True;
        }

        // press q to quit
        else if (key == 'q')            exit(0);

        // toggle rotating
        else if (key == 'r') {
            if (rotate == True) {
                rotate    = False;
                translate = True;
            }
            else {
                rotate    = True;
                translate = False;
            }
        }

        // toggle translating
        else if (key == 't') {
            if (translate == True) {
                rotate    = True;
                translate = False;
            }
            else {
                rotate    = False;
                translate = True;
            }
        }

        // make translation matrix and translate
        else if (translate == True) {
            if (key == 'x') {
                if      (reverse == True)   M3d_copy_mat(transform, x_trans_back);
                else if (reverse == False)  M3d_copy_mat(transform, x_trans_forw);
            }
            else if (key == 'y') {
                if      (reverse == True)   M3d_copy_mat(transform, y_trans_back);
                else if (reverse == False)  M3d_copy_mat(transform, y_trans_forw);
            }
            else if (key == 'z') {
                if      (reverse == True)   M3d_copy_mat(transform, z_trans_back);
                else if (reverse == False)  M3d_copy_mat(transform, z_trans_forw);
            }
        }

        // make rotation matrix
        else if (rotate == True) {
            // translate_center * rotate * translate_back
            M3d_make_translation(transform, -x[file_num][num_points[file_num]], -y[file_num][num_points[file_num]], -z[file_num][num_points[file_num]]);

            if (key == 'x') {
                if      (reverse == True)   M3d_mat_mult(transform, x_rot_back, transform);
                else if (reverse == False)  M3d_mat_mult(transform, x_rot_forw, transform);
            }
            else if (key == 'y') {
                if      (reverse == True)   M3d_mat_mult(transform, y_rot_back, transform);
                else if (reverse == False)  M3d_mat_mult(transform, y_rot_forw, transform);
            }
            else if (key == 'z') {
                if      (reverse == True)   M3d_mat_mult(transform, z_rot_back, transform);
                else if (reverse == False)  M3d_mat_mult(transform, z_rot_forw, transform);
            }

            // this should work as if it didnt fit other conditions it will continue anyways?
            else {
                continue;
            }

            // now translate back to original location?
            double translate_back[4][4];
            M3d_make_translation(translate_back, x[file_num][num_points[file_num]], y[file_num][num_points[file_num]], z[file_num][num_points[file_num]]);
            M3d_mat_mult(transform, translate_back, transform);
        }

        else {
            continue;
        }

        // multiply matrix by points
        M3d_mat_mult_points(x[file_num], y[file_num], z[file_num], transform, x[file_num], y[file_num], z[file_num], num_points[file_num] + 1);

        G_rgb(0, 0, 0);
	    G_fill_rectangle(0, 0, window_size, window_size);

        printf("Reverse %d\nRotate %d\nTranslate %d\n", reverse, rotate, translate);
        printf("Displaying %d\n\n", file_num);
        G_rgb(0, 1, 0);
        draw_polygon(file_num);
    }
}
