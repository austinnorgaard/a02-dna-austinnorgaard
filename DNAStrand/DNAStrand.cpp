//***************************************************************************
// File name:		DNAStrand.cpp
// Author:			Austin Norgaard
// Date:				02/07/2022
// Class:				CS 250
// Assignment:  01 DNAStrand Class
// Purpose:			Define the interface for a strand of DNA. A strand of DNA
//							contains a unique unknown number of chemical bases
//							adenine (A), guanine (G), cytosine (C), and thymine (T).
//							The ID will be a simplified FASTA string.
//***************************************************************************

#include <iostream>
#include <fstream>
#include "DNAStrand.h"

using namespace std;

/****************************************************************************
Function:			DNAStrand Default Constructor

Description:	The default constructor for the DNAStrand object

Parameters:		rcID				 - The ID used for the initialization of a 
														 DNAStrand, default = ""
							rcDNAStrand  - The Bases used for the initialization of a 
														 DNAStrand, default = ""

Returned:			none
****************************************************************************/
DNAStrand::DNAStrand (const string &rcID, const string &rcDNAStrand)
{
	mID = rcID;
	mBases = rcDNAStrand;
}

/****************************************************************************
Function:			DNAStrand Copy Constructor

Description:	A copy constructor for the DNAStrand object

Parameters:		rcDNAStrand  - The DNAStrand object in which the values are
														 being copied from

Returned:			none
****************************************************************************/
DNAStrand::DNAStrand	(const DNAStrand &rcDNAStrand)
{
	mID = rcDNAStrand.mID;
	mBases = rcDNAStrand.mBases;
}

/****************************************************************************
Function:			getBase

Description:	Return the Base of a DNAStrand at a particular location 

Parameters:		whichBase  - The position of the char being returned

Returned:			base       - The base of the DNAStrand at the position passed
													 in
****************************************************************************/
char DNAStrand::getBase (int whichBase) const
{
	char base = mBases[whichBase];

	return base;
}

/****************************************************************************
Function:			read

Description:	Read data from an input stream into the members of the
							DNAStrand as an ID then the Bases

Parameters:		rcInStream				 - The input stream being read from

Returned:			LEGAL_BASE_FLAG	   - Returned if the reading is done
																	 correctly and there are no errors
							ILLEGAL_BASE_FLAG  - Returned if the reading does not finish
																	 and there are illegal bases read into
																	 mBases
****************************************************************************/
enum DNAStrand::FLAG DNAStrand::read (istream &rcInStream)
{
	const int LENGTH = 3;
	// Verify the file contains data for mID and mBases for at least 1 strand
	if (rcInStream >> mID >> mBases)
	{
		for (int i = 0; i < mBases.length (); i++)
		{
			// Verify all bases are legal
			if (mBases[i] != ADENINE && mBases[i] != CYTOSINE &&
					mBases[i] != GUANINE && mBases[i] != THYMINE)
			{
				ILLEGAL_BASE;
			}
			if (mBases == "" || mID == "")
			{
				ILLEGAL_BASE;
			}
			if (mBases.size () != LENGTH)
			{
				ILLEGAL_BASE;
			}
		}
	}
	LEGAL;
}

/****************************************************************************
Function:			write

Description:	Write the members of the DNAStrand object to an output stream

Parameters:		rcOutStream  - The output stream being written to

Returned:			none
****************************************************************************/
void DNAStrand::write (ostream &rcOutStream) const
{
	rcOutStream << mID << " ";

	writeBases (rcOutStream);

	rcOutStream << endl;
}

/****************************************************************************
Function:			writeBases

Description:	Write the bases of the DNAStrand object only, to an output 
							stream

Parameters:		rcOutStream  - The output stream being written to

Returned:			none
****************************************************************************/
void DNAStrand::writeBases (ostream &rcOutStream) const
{
	rcOutStream << mBases;
}

/****************************************************************************
Function:			reverse

Description:	Return the reversed bases of a DNAStrand

Parameters:		rcID     - The ID of the returned DNAStrand

Returned:			cStrand  - The reversed DNAStrand
****************************************************************************/
DNAStrand DNAStrand::reverse (const string &rcID) const
{
	DNAStrand cStrand;
	// Static Cast size_t to int to prevent loss of data
	int length = static_cast<int> (mBases.size ());
	// Initialize mBases to prevent possible out of range error
	cStrand.mBases = mBases;

	for (int i = length - 1, k = 0; i > -1; i--, k++)
	{
		cStrand.mBases[k] = mBases[i];
	}

	cStrand.mID = rcID;

	return cStrand;
}

/****************************************************************************
Function:			complement

Description:	Return the complemented bases of a DNAStrand

Parameters:		rcID     - The ID of the returned DNAStrand

Returned:			cStrand  - The complemented DNAStrand
****************************************************************************/
DNAStrand DNAStrand::complement (const string &rcID) const
{
	DNAStrand cStrand;
	// Initialize mBases to prevent possible out of range error
	cStrand.mBases = mBases;

	for (int i = 0; i < mBases.size (); i++)
	{
		cStrand.mBases[i] = dnaBaseComplement (mBases[i]);
	}
	
	cStrand.mID = rcID;

	return cStrand;
}

/****************************************************************************
Function:			gcContent

Description:	Return the GC Content of a DNAStrand

Parameters:		none

Returned:			gcContent  - The ratio of G and C within a DNAStrands bases
													 compared to the total number of bases
****************************************************************************/
double DNAStrand::gcContent	() const
{
	double gcContent = 0.0;

	for (int i = 0; i < mBases.size	(); i++)
	{
		if (mBases[i] == CYTOSINE || mBases[i] == GUANINE)
		{
			gcContent += 1.0;
		}
	}

	gcContent /= mBases.size	();

	return gcContent;
}

/****************************************************************************
Function:			hamming

Description:	Return the Hamming Distance of two DNAStrands

Parameters:		none

Returned:			hammingDistance  - The difference of bases between two
																 DNAStrands compared to the total length
																 of a DNAStrand
****************************************************************************/
double DNAStrand::hamming	(const DNAStrand &rcDNAStrand) const
{
	double hammingDistance = 0.0;

	for (int i = 0; i < mBases.size (); i++)
	{
		if (mBases[i] != rcDNAStrand.mBases[i])
		{
			hammingDistance += 1.0;
		}
	}

	hammingDistance /= mBases.size ();

	// Verify size of mBases of both DNA
	if (mBases.size () != rcDNAStrand.mBases.size ())
	{
		hammingDistance = ILLEGAL_EOF;
	}

	return hammingDistance;
}

/****************************************************************************
Function:			size

Description:	Return the size of a DNAStrand

Parameters:		none

Returned:			length  - The total number of bases within a DNAStrand
****************************************************************************/
unsigned int DNAStrand::size () const
{
	unsigned int length = 0;

	while (length < mBases.size ()) 
	{
		length++;
	}

	return length;
}

/****************************************************************************
Function:			isEqual

Description:	Return whether the IDs of two DNAStrands are the same.

Parameters:		rcDNAStrand  - The DNAStrand being compared to

Returned:			isEqual      - Either True or False depending on if the two
														 strands IDs are equal
****************************************************************************/
bool DNAStrand::isEqual (const DNAStrand &rcDNAStrand) const
{
	return rcDNAStrand.mID == mID;
}

/****************************************************************************
Function:			dnaBaseComplement

Description:	The complement of each base

Parameters:		base  - The base being complemented

Returned:			base  - The complemented base
****************************************************************************/
char DNAStrand::dnaBaseComplement (char base) const
{
	// A = T, T = A, C = G, G = C
	if (base == ADENINE)
	{
		base = THYMINE;
	}
	else if (base == THYMINE)
	{
		base = ADENINE;
	}
	else if (base == CYTOSINE)
	{
		base = GUANINE;
	}
	else if (base == GUANINE)
	{
		base = CYTOSINE;
	}

	return base;
}