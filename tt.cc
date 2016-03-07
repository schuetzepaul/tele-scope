
// Daniel Pitzl, Feb 2016
// test file tokens

// g++ -o tt tt.cc

#include <iostream> // cout
#include <iomanip> // setw
#include <fstream>
#include <sstream>

using namespace std;

int main()
{
  ifstream tfile( "tt.dat" );

  if( tfile.bad() || !tfile.is_open() ) {
    cout << "Error opening tt.dat" << endl;
    return 1;
  }

  string hash( "#" );

  int nlines = 0;

  while( ! tfile.eof() ) {

    string line;
    getline( tfile, line );
    cout << setw(3) << nlines << ": " << line << endl;
    ++nlines;
    if( line.empty() ) {
      cout << "empty" << endl;
      continue;
    }
    stringstream tokenizer( line ); // break into whitespaced blocks
    string tag;
    tokenizer >> tag; // leading whitespace is dropped
    cout << "tag: " << tag << endl;
    if( tag.substr(0,1) == hash ) {
      cout << "ignore" << endl;
      continue;
    }
    int val;
    tokenizer >> val;
    cout << "val: " << val << endl;

  } // eof

  cout << endl << "eof after " << nlines << " lines" << endl << endl;

  return 0;
  }
