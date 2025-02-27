#include "rtwtypes.h"

void change_gpio_value(int gpioPin, int newValue);
void export_gpio(int gpioPin);
void unexport_gpio(int gpioPin);
void write_gpio_value(int gpioPin, int value);
void set_pin_direction(int gpioPin, int direction);
