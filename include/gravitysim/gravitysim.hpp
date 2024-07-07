#pragma once

struct Spaceobject {
    // units in m
    double x;
    double y;

    // units in m/s
    double speedX = 0;
    double speedY = 0;

    // radius in m
    double radius;
    // mass in kg (guesstimated as sphere with density of 2g*cm^-1)
    double mass;

    Spaceobject(double _x, double _y) : x(_x), y(_y) {
        radius = 500;
        mass = (1.33 * 3.1415 * radius * radius * radius) * 2;
    }
    Spaceobject(double _x, double _y, double _radius) : x(_x), y(_y), radius(_radius) {
        mass = (1.33 * 3.1415 * radius * radius * radius) * 2;
    }
};