// Light modeling

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <FPT.h>
#include "M3d_matrix_toolsS.c"

double light[3];
double ambient;
double diffuse_max;

int window_size = 800;
int num_points[10000], num_polys[10];
int psize[50][1000];
int con[100][1000][1000];
double x[10][50000], y[10][50000], z[10][50000];

typedef
struct {
    int poly_num; // polygon
    int obj_num;  // object poly is in
    double distance; // distance from origin
}
POLYGON;
double dot_product(double vector1[3], double vector2[3]) {
    return ((vector1[0] * vector2[0]) + (vector1[1] * vector2[1]) + (vector1[2] * vector2[2]));
}

double cross_product(double result[3], double vector1[3], double vector2[3]) {
    result[0] =   (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
    result[1] = -((vector1[0] * vector2[2]) - (vector1[2] * vector2[0]));
    result[2] =   (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);
}

int make_vector(double vector[3], double x1, double y1, double z1, double x2, double y2, double z2){
    vector[0] = x2 - x1;
    vector[1] = y2 - y1;
    vector[2] = z2 - z1;
}

int make_perp_vector(double perp_vector[3],
		     double x1, double y1, double z1,
		     double x2, double y2, double z2,
		     double x3, double y3, double z3) {
    double vector1[3], vector2[3];
    make_vector(vector1, x1, y1, z1, x2, y2, z2);
    make_vector(vector2, x1, y1, z1, x3, y3, z3);

    cross_product(perp_vector, vector1, vector2);
}

int compare(const void *p, const void *q) {
    POLYGON *a, *b;
    a = (POLYGON*)p;
    b = (POLYGON*)q;

    if      (((*a).distance) < ((*b).distance)) return  1;
    else if (((*a).distance) > ((*b).distance)) return -1;
    else                                        return  0;
}

int find_poly_center(double *center, double *x, double *y, double *z, int length) {
    double x_total = 0, y_total = 0, z_total = 0;

    // gets the x, y, z values of each point that is in the polygon
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

    // need a count, as number of points in entire polygon isn't certain
    int count = 0;

    // get the x, y, z values of each point that is in each polygon
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

double find_length(double vector[3]) {
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}

double find_intensity(double x_points[5000], double y_points[5000], double z_points[5000], int length) {
    double center[3];
    find_poly_center(center, x_points, y_points, z_points, length);

    // make view vector, light vector, perp vector
    // I know light is spelled wrong, I did it so it would be the same length as view and perp
    double view_vector[3], lite_vector[3], perp_vector[3];
    make_vector(view_vector,    x_points[0], y_points[0], z_points[0],
                                0, 0, 0);
    make_vector(lite_vector,    x_points[0], y_points[0], z_points[0],
                                light[0],  light[1],  light[2]);
    make_perp_vector(perp_vector,   x_points[0], y_points[0], z_points[0],
                                    x_points[1], y_points[1], z_points[1],
                                    x_points[2], y_points[2], z_points[2]);

    // make these vectors into unit vectors!
    double view_length = find_length(view_vector);
    view_vector[0] /= view_length;
    view_vector[1] /= view_length;
    view_vector[2] /= view_length;

    double lite_length = find_length(lite_vector);
    lite_vector[0] /= lite_length;
    lite_vector[1] /= lite_length;
    lite_vector[2] /= lite_length;

    double perp_length = find_length(perp_vector);
    perp_vector[0] /= perp_length;
    perp_vector[1] /= perp_length;
    perp_vector[2] /= perp_length;

    // If light and eye are on opposite sides, intensity = ambinet
    //    double dot_lite_view = dot_product(lite_vector, view_vector);
    double dot_lite_perp = dot_product(lite_vector, perp_vector);
    double dot_view_perp = dot_product(view_vector, perp_vector);

    // if light is on the other side of the polygon, just use ambient light
    if (dot_lite_perp * dot_view_perp < 0) {
        return ambient;
    }

    // if light and eye are on opposite side of normal vector, reverse perp vector
    else if (dot_lite_perp < 0 && dot_view_perp < 0) {
        perp_vector[0] *= -1;
        perp_vector[1] *= -1;
        perp_vector[2] *= -1;

	dot_lite_perp *= -1 ;
	dot_view_perp *= -1 ;	
	
    }

    // find reflection vector, Ru = -Lu + 2(Nu . Lu) . Nu
    double refl_vector[3];
    refl_vector[0] = -lite_vector[0] + 2 * (dot_lite_perp * perp_vector[0]);
    refl_vector[1] = -lite_vector[1] + 2 * (dot_lite_perp * perp_vector[1]);
    refl_vector[2] = -lite_vector[2] + 2 * (dot_lite_perp * perp_vector[2]);

    //printf("Ambient %lf\n", ambient);
    //TODO issue with perp and lite vectors here. Is this because of the values? Math seems to check out

    //
    //printf("Dot %lf\n", (dot_product(perp_vector, lite_vector)));
    //printf("Perp %lf, %lf, %lf, lite %lf, %lf, %lf\n", perp_vector[0], perp_vector[1], perp_vector[2], lite_vector[0], lite_vector[1], lite_vector[2]);
    //
    // printf("View/refl %lf\n", viewrefl);

    double intensity = ambient + diffuse_max * dot_lite_perp + (1 - ambient - diffuse_max) * pow(dot_product(view_vector, refl_vector),50);

    if (intensity > 1 || intensity < 0) {
        printf("Intensity %lf\n", intensity);
        double perplite = (diffuse_max * dot_product(perp_vector, lite_vector));
        double viewrefl = ((1 - ambient - diffuse_max) * dot_product(view_vector, refl_vector));
        printf("Perp/lite %lf\n", perplite);
        printf("View/refl %lf\n", viewrefl);

        printf("\n");
    }

    return intensity;
}

int set_rgb(double x_points[5000], double y_points[5000], double z_points[5000], int length) {
    double intensity = find_intensity(x_points, y_points, z_points, length);

    double red = 0.75;
    double gre = 0;
    double blu = 0.75;

    double constants = ambient + diffuse_max;

    if (intensity < (ambient + diffuse_max)) {
        red *= intensity / constants;
        gre *= intensity / constants;
        blu *= intensity / constants;
    }
    // is this the right equation?
    else if (intensity > (ambient + diffuse_max)) {
        red *= ((intensity - constants) / (1 - constants)) + 1;
        gre *= ((intensity - constants) / (1 - constants)) + 1;
        blu *= ((intensity - constants) / (1 - constants)) + 1;
    }

    printf("RGB %lf %lf %lf\n\n", red, gre, blu);
    G_rgb(red, gre, blu);
    //G_rgb(intensity,intensity,intensity);
}

void draw_polygon(int num_objects) {
    int total_poly_num = 0;

    // calculate total number of polygons being drawn
    for (int obj_num = 0; obj_num < num_objects; obj_num++) {
        total_poly_num += num_polys[obj_num];
    }

    // creates array of polygons
    POLYGON polygons[total_poly_num];

    int count = 0;

    // calculates each polygon's distance from the viewer
    for (int obj_num = 0; obj_num < num_objects; obj_num++) {
        // calculate each polygon in each object
        for (int poly_num = 0; poly_num < num_polys[obj_num]; poly_num++) {
            double x_points[1000], y_points[1000], z_points[1000];

            // get the x, y values of each point that is in each polygon
            for (int point = 0; point < psize[obj_num][poly_num]; point++) {
                x_points[point] = x[obj_num][con[obj_num][poly_num][point]];
                y_points[point] = y[obj_num][con[obj_num][poly_num][point]];
                z_points[point] = z[obj_num][con[obj_num][poly_num][point]];
            }

            // do math to calculate 3d aspect of this shape
            double x_bar, y_bar, x_barbar, y_barbar, h;
            h = tan(60);

            for (int point = 0; point < psize[obj_num][poly_num]; point++) {
                x_bar = x_points[point] / z_points[point];
                y_bar = y_points[point] / z_points[point];

                x_barbar = (window_size / 2) * (x_bar / h) + (window_size / 2);
                y_barbar = (window_size / 2) * (y_bar / h) + (window_size / 2);

                x_points[point] = x_barbar;
                y_points[point] = y_barbar;
            }

            // find center of polygon, make distance the z value of the center
            double center[3];
            find_poly_center(center, x_points, y_points, z_points, psize[obj_num][poly_num]);

            polygons[count].obj_num  = obj_num;
            polygons[count].poly_num = poly_num;
            polygons[count].distance = center[2];

            count += 1;
        }
    }

    // sorts polygon with library quicksort
    qsort(polygons, total_poly_num, sizeof(POLYGON), compare);

    // draws polygons farthest, first
    for (int poly = 0; poly < total_poly_num; poly++) {
        double x_points[1000], y_points[1000], z_points[1000];

        // sets polygon values
	int obj_num, poly_num;
	double distance;
        obj_num  = polygons[poly].obj_num;
        poly_num = polygons[poly].poly_num;
        distance = polygons[poly].distance;

        //printf("Obj_num %d, poly_num %d\n", obj_num, poly_num);

        int size = psize[obj_num][poly_num];

	    //printf("Before setting x, y, z points\n");

        // get the x, y values of each point that is in each polygon
        for (int point = 0; point < size; point++) {
            x_points[point] = x[obj_num][con[obj_num][poly_num][point]];
            y_points[point] = y[obj_num][con[obj_num][poly_num][point]];
            z_points[point] = z[obj_num][con[obj_num][poly_num][point]];
        }


        set_rgb(x_points, y_points, z_points, psize[obj_num][poly_num]);
	
        // do math to calculate 3d aspect of this shape
        double x_bar, y_bar, x_barbar, y_barbar, h;
        h = tan(45*M_PI/180);

        for (int point = 0; point < psize[obj_num][poly_num]; point++) {
            x_bar = x_points[point] / z_points[point];
            y_bar = y_points[point] / z_points[point];

            x_barbar = (window_size / 2) * (x_bar / h) + (window_size / 2);
            y_barbar = (window_size / 2) * (y_bar / h) + (window_size / 2);

            x_points[point] = x_barbar;
            y_points[point] = y_barbar;
        }

	

        // print polygon

        //G_rgb(0, 1, 0);
        G_fill_polygon(x_points, y_points, psize[obj_num][poly_num]);

        //G_polygon(x_points, y_points, psize[obj_num][poly_num]);
    }
}

int main(int argc, char **argv) {
  light[0] = 100 ;
  light[1] = 200 ;
  light[2] = 0 ;
    ambient = 0.2;
    diffuse_max = 0.5;

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
        // storing center in x, y, z, so that operations on the figure will also affect the center
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

    draw_polygon(num_objects);

    // loop for getting keys and drawing polygons
    while (1) {
        double transform[4][4];
        M3d_make_identity(transform);

        int key = G_wait_key();
        int key_num = key - 49;

        printf("Key_num is %d\n", key_num);

        // display new object
        if (key_num < num_objects)  file_num = key_num;

        // reverse the transform/rotation
        else if (key == 'c') {
            if (reverse == True)    reverse = False;
            else                    reverse = True;
        }

        // press q to quit
        else if (key == 'q')        exit(0);

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

            // this should work as if it didnt fit other conditions it will continue anyways
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

        // clear screen
        G_rgb(0, 0, 0);
	    G_fill_rectangle(0, 0, window_size, window_size);

        printf("Reverse %d\nRotate %d\nTranslate %d\n", reverse, rotate, translate);
        printf("Displaying %d\n\n", file_num);
        G_rgb(0, 1, 0);

        // draws objects
        draw_polygon(num_objects);
        //draw_polygon(2);
    }
}
