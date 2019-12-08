#include <FPT.h>

// compares values for quicksort
// @param a is the first value
// @param b is the second value
int compare(const void * a, const void * b) {
	return ( *(double*)a - *(double*)b );
}

// has to fill in all of the points in a polygon
// @param x[] is the x coordinates
// @param y[] is the y coordinates
// @param n is the number of points
int G_fill_my_polygon(double x[], double y[], int n) {
	// iterate through each x value
	// if a line intersects that x value, start coloring in pixels
	// if another line intersects that x value, stop coloring in pixels		
	for (int y_value = 0; y_value <= 400; y_value++) {
		double intersections[n];
		double previous_y = y[n - 1];
		int count = 0;

		// iterate through all vertices
		// looks at vertices it's connected to and determine where it intercepts this y_value
		for (int vertex = 0; vertex < n; vertex += 1) {
			// sets the x and y values for v1
			double v1x = x[vertex];
			double v1y = y[vertex];

			// sets v2 values to next point, unless vertex is the last in array
			double v2x, v2y;
			if (vertex < (n - 1)) {
				v2x = x[vertex + 1];
				v2y = y[vertex + 1];
			}
			else {
				v2x = x[0];
				v2y = y[0];
			}

			// find which y value is bigger
			double bigger_y, smaller_y;
			if (v1y > v2y) {
				bigger_y = v1y;
				smaller_y = v2y;
			}
			else {
				bigger_y = v2y;
				smaller_y = v1y;
			}

			// find if current y value is out of bounds of this line
			if (bigger_y < y_value || smaller_y > y_value ) {
				continue;
			}

			// find equation for line
			double slope = (v2y - v1y) / (v2x - v1x);
			double y_intercept = (slope * v1x) - v1y;
			
			// if y = y_intercept, add each x_value
			if (slope == 0) {
				continue;
			}

			double x_value = (y_intercept + y_value) / slope;

			// if x = x_intercept, make x_value x_intercept
			if (v1x == v2x) {
				x_value = v1x;
			}

			// add x value
			intersections[count] = x_value;
			count += 1;
			
			// check for non tips and double those entries
			if (v1y == y_value) {
				if ((v1y > v2y && v1y < previous_y) || (v1y < v2y && v1y > previous_y)) {
					intersections[count] = x_value;
					count += 1;
				}
			}

			previous_y = v1y;
		}	
		// stdlib qsort
		qsort(intersections, count, sizeof(double), compare);

		// draw line from each vertex to the next, skipping every other
		for (int vertex = 0; vertex < count; vertex += 1) {
			if (vertex % 2 == 0) { 
				G_line(intersections[vertex], y_value, intersections[vertex + 1], y_value);
			}
		}
	}
}


// allows user to click multiple points, and then when the user clicks the yellow rectangle the points are filled in
int main() {
	double x[100], y[100];

	G_init_graphics(400, 400);
	
	// creates yellow rectangle for user to click to end their polygon drawing	
	G_rgb(1, 1, 0);
	G_fill_rectangle(0, 0, 50, 50);

	G_rgb(0, 0, 0);


	// loop for user to click the points for their polygon
	int count = 0;

	while(1) {
		double p[2];
		G_wait_click(p);

		// breaks if user clicks where the yellow rectangle is
		if ((p[0] <= 50) && (p[1] <= 50)) {
			break;
		}
		else {
			G_point(p[0], p[1]);

			x[count] = p[0];
			y[count] = p[1];

			if (count > 0) {
				G_line(x[count - 1], y[count - 1], x[count], y[count]);
			}
			count += 1;
		}
	}

	G_fill_my_polygon(x, y, count);

	G_wait_key();
}
