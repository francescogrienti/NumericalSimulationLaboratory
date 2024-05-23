//
// Created by francesco on 23/05/24.
//

#ifndef NUMERICALSIMULATIONLABORATORY_CITY_H
#define NUMERICALSIMULATIONLABORATORY_CITY_H


class City {
private:
    int city_label;
    double angular_coordinate;
public:
    // Default constructor
    City();

    // Destructor
    ~City();

    //Method for getting the label of the city
    int getLabel();

    //Method for getting the angular coordinate of the city
    double getCoordinate();

    //Method for setting the label of the city
    void setLabel(int label);

    //Method for setting the angular coordinate of the city
    void setCoordinate(double coordinate);

};


#endif //NUMERICALSIMULATIONLABORATORY_CITY_H
