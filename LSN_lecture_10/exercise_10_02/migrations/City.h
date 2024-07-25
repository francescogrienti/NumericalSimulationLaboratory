//
// Created by francesco on 23/05/24.
//

#ifndef NUMERICALSIMULATIONLABORATORY_CITY_H
#define NUMERICALSIMULATIONLABORATORY_CITY_H

#include <vector>

class City {
private:
    int city_label;
    double angular_coordinate; //Circle
    double x, y; //Square
public:
    // Default constructor
    City();

    // Destructor
    ~City();

    //Method for getting the label of the city
    int getLabel();

    //Method for getting the angular coordinate of the city
    double getCoordinate();

    //Method for getting the coordinates of the city in the square
    std::vector<double> getSquareCoordinates();

    //Method for setting the label of the city
    void setLabel(int label);

    //Method for setting the angular coordinate of the city
    void setCoordinate(double coordinate);

    //Method for setting the angular coordinate of the city
    void setSquareCoordinates(double x_1, double x_2);
};


#endif //NUMERICALSIMULATIONLABORATORY_CITY_H
