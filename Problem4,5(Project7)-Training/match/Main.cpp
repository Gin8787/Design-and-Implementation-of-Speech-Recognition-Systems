#include "gaussian-DTW.h"
#include "ktrain.h"

string Template, Input;

int main() {
	int i;
	for (i = 0; i < 10; i++){
		segKMeansProcess();
		assigndigits();
		clearall();
	}
	return 0;
}