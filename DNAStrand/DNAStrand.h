#pragma once
//********************************************************************* 
// File name:		DNAStrand.h
// Author:			CS, Pacific University
// Date:				
// Class:				CS 250
// Assignment:  01 DNA Class
// Purpose:			Declare the interface for a strand of DNA. A strand of
//							DNA contains a unique and an unknow number of chemical
//							bases adenine (A), guanine (G), cytosine (C), and 
//							thymine (T). The ID will be a simplified FASTA string.
//*********************************************************************


#include <string>
#include <iostream>

using namespace std;

class DNAStrand {

public:
	enum FLAG : int {
		ILLEGAL_BASE = -2, ILLEGAL_EOF = -1, LEGAL = 0, LEGAL_EOF = 1
	};

	enum BASE : char {
		ADENINE = 'A', CYTOSINE = 'C', GUANINE = 'G',
		THYMINE = 'T', NO_BASE = ' '
	};

	DNAStrand (const string &rcID = "", const string &rcBases = "");
	DNAStrand (const DNAStrand &rcDNAStrand);

	char getBase (int whichBase) const;

	FLAG read (istream &rcInStream);
	void write (ostream &rcOutStream) const;
	void writeBases (ostream &rcOutStream) const;

	DNAStrand reverse (const string &rcID) const;
	DNAStrand complement (const string &rcID) const;
	double gcContent () const;
	double hamming (const DNAStrand &rcDNAStrand) const;

	unsigned int size () const;
	bool isEqual (const DNAStrand &rcDNAStrand) const;

private:
	string mID;
	string mBases;

	char dnaBaseComplement (char base) const;
};