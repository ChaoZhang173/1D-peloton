#ifndef __INITIALIZER_H__
#define __INITIALIZER_H__

using namespace std;
class Global_Data;

/*
This file is for initializer class. 
It is used to initialize the parameters and pass them to Global_Data class.
In later version may use a input file to initialize, but currently only set to default values.

*/

class Initializer {
friend class Global_Data; // Global_Data has access to private and proteced members
public:
// constructor
Initializer();

private:
// currently not using
void readInputfile(const std::string& inputfileName);
void setParams();

}
