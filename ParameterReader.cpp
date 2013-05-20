/***********************************************************************
See ParameterReader.h for brief description and change log.
09-09-2011 Zhi Qiu
***********************************************************************/
#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "arsenal.h"
#include "ParameterReader.h"

using namespace std;

//----------------------------------------------------------------------
ParameterReader::ParameterReader()
{
  names = new vector<string>;
  values = new vector<double>;
  paths = new vector<string>;
}


//----------------------------------------------------------------------
ParameterReader::~ParameterReader()
{
  delete names; // in principle garbage collection for vector should be automatic, just to be safe
  delete values;
  delete paths;
}


//----------------------------------------------------------------------
string ParameterReader::removeComments(string str, string commentSymbol)
/*
  Remove comments from a string "str". Comments are all characters after the string "commentSymbol".
*/
{
  return str.substr(0, str.find(commentSymbol));
}


//----------------------------------------------------------------------
void ParameterReader::phraseEquationWithoutComments(string equation)
/*
  Phrase an equation like "x=1", and store the result into "names" and "values". The equation is first separated according to the equal sign, then the left and right hand side will be trimmed, after that the right hand side will be converted to double type number.
*/
{
  if (trim(equation).compare("")==0) return;
  size_t symbolPos = equation.find('=');
  if (symbolPos==string::npos)
  {
    cout << "ParameterReader::phraseEquationWithoutComments error: \"=\" symbol not found in equation assignment " << equation << endl;
    exit(-1);
  }
  string LHS (equation.begin(), equation.begin()+symbolPos);
  string RHS (equation.begin()+symbolPos+1, equation.end());
  if(trim(LHS)=="path") setPath(LHS, trim(RHS)); 
  else setVal(LHS, stringToDouble(trim(RHS))); 
}


//----------------------------------------------------------------------
long ParameterReader::find(string name)
/*
  Check if the parameter with "name" already exists in the internal "names" list. If yes, it returns its
*/
{
  for (long ii=0; ii<names->size(); ii++)
    if ((*names)[ii].compare(toLower(trim(name)))==0) return ii;
  return -1;
}


//----------------------------------------------------------------------
void ParameterReader::phraseOneLine(string str, string commentSymbol)
/*
  Interpret a string like " x  = 1.1  #bla " to get the associated parameter name and value information, and put them into the internal variables "names" and "values".
*/
{
  if (trim(str).compare("")==0) return;
  phraseEquationWithoutComments(removeComments(str, commentSymbol));
}


//----------------------------------------------------------------------
void ParameterReader::readFromFile(string filename, string commentSymbol)
/*
  Read all lines in a file as parameter assignment list. Each line is processed by the phraseOneLine function.
*/
{
  ifstream parameterFile(filename.c_str());
  if (!parameterFile)
  {
    cout << "ParameterReader::readFromFile error: file " << filename << " does not exist." << endl;
    exit(-1);
  }
  char buffer[9999];
  while (!parameterFile.eof())
  {
    parameterFile.getline(buffer, 9999);
    phraseOneLine(buffer);
  }
  parameterFile.close();
}


//----------------------------------------------------------------------
void ParameterReader::readFromArguments(long argc, char * argv[], string commentSymbol, long start_from)
/*
  Read all strings in argv[]. Each string is processed by the phraseOneLine function.
*/
{
  for (long ii=start_from; ii<argc; ii++) phraseOneLine(argv[ii], commentSymbol);
}


//----------------------------------------------------------------------
bool ParameterReader::exist(string name)
/*
  Return true if parameter with "name" is registered.
*/
{
  return find(name)==-1 ? false: true;
}


//----------------------------------------------------------------------
void ParameterReader::setVal(string name, double value)
/*
  Set the parameter with "name" to "value". It is appended to the internal "names" and "values" vector if "name" does not exist; otherwise it is rewitten.
*/
{
  long idx = find(name);
  if (idx==-1)
  {
    names->push_back(toLower(trim(name))); 
    values->push_back(value);
    paths->push_back("null");
  }
  else
  {
    (*names)[idx]=toLower(trim(name)); (*values)[idx]=value; (*paths)[idx]="null";
  }
}

void ParameterReader::setPath(string name, string path)
/*
  Set the parameter with "name" to "path". It is appended to the internal "names" and "paths" vector if "name" does not exist; otherwise it i  s written.
*/
{
  long idx = find(name);
  if (idx==-1)
  {
    names->push_back(toLower(trim(name)));
    values->push_back(stringToDouble(path));
    paths->push_back(path);
  }
  else
  {
    (*names)[idx]=toLower(trim(name)); (*values)[idx]=stringToDouble(path); (*paths)[idx]=path;
  }
}

//----------------------------------------------------------------------
double ParameterReader::getVal(string name)
/*
  Get the value for the parameter with "name".
*/
{
  long idx = find(name);
  if (idx!=-1)
    return (*values)[idx];
  else
  {
    cout << "ParameterReader::getVal error: parameter with name " << name << " not found." << endl;
    exit(-1);
  }
}

string ParameterReader::getPath(string name)
/*
  Get the path for the parameter with "path".
*/
{
  long idx = find(name);
  if (idx!=-1)
    return (*paths)[idx];
  else
  {
    cout << "ParameterReader::getPath error: parameter with name " << name << " not found." << endl;
    exit(-1);
  }
}


//----------------------------------------------------------------------
void ParameterReader::echo()
/*
  Print out all stored parameters to screen.
*/
{
  if (names->size()==0) return;
  for (long ii=0; ii<names->size(); ii++){
	if((*names)[ii]=="path") cout << (*names)[ii] << "=" << (*paths)[ii] << "  "; 
	else cout << (*names)[ii] << "=" << (*values)[ii] << "  ";
	}
  cout << endl;
}
