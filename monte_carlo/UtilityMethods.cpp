/*
 * UtilityMethods.c++
 *
 */
#include "UtilityMethods.h"

double computeDistanceSquare(Vector3d position1,Vector3d position2){
	   //  distance between two vectors

		double distance;
		Vector3d tempVector;

		tempVector.x = position2.x - position1.x ;
		tempVector.y = position2.y - position1.y ;
		tempVector.z = position2.z - position1.z ;

		distance = tempVector.x*tempVector.x + tempVector.y*tempVector.y +
				   tempVector.z*tempVector.z ;

		return distance;
}

KahanSum::KahanSum(void){
	sum = 0.0;
	correction = 0.0;
}

void KahanSum::add(double v){
	 double y = v - correction;
	 double t = sum + y;
	 // (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
	 correction = (t - sum) - y;
	 sum = t; //Next time around, the lost low part c will be added to
}

void KahanSum::add(vector<double> v){
	for (auto d : v)
			add(d);
}
double KahanSum::getSum(void){
	return sum;
}
