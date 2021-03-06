/* Rooter: stores MOL output from stdin in a ROOT TTree of 64-bit values
 * 
 * Syntax: cat input.dat | rooter output.root
 */
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <iostream>
#include <stdint.h>

#define GZ 0x5a47 // FEROL "magic number"
#define START (1L<<39) // initial block: bit 31 --> 39 (endian flip), as well as zeroes in 42..32
#define FEDid 931 //Define FEDid to add into the branch name
#define SaveVecBranch 1 //Switch to choose between storing vector or array 
#define SaveHeader 0//Switch to choose if we want to store the header word into the root file

int DecodeNDataInBlock(uint64_t buf){
  int nData = (buf>>56);
  nData += (buf>>40) & 0x300;
  return nData;
}

int main(int argc, char* argv[]) {
  gROOT->ProcessLine("#include <vector>; #pragma link C++ class vector<uint64_t>+;"); //Load dictionary for vector<uint64_t>
  char* rootfile;
  if (argc == 1) strcpy(rootfile, "mol.root");
  else rootfile = argv[1];
  TFile a(rootfile, "RECREATE");
  TTree tree("moltree","");
  
  const int MAX_WORDS = 100000;

  uint64_t blob[MAX_WORDS];
  std::vector<uint64_t> vec(MAX_WORDS);
  vec.clear();
  int len = 0;
  if(SaveVecBranch) tree.Branch(Form("vec%d",FEDid), &vec);
  else{
    tree.Branch("nWord64",&len,"nWord64/I");
    tree.Branch("words",blob,"words[nWord64]/l");
  }
  uint64_t buf;
  int n = 0; // fragments
  int iWordInBlock = 0;
  int nWordsInBlock = 0;
  int success = fread(&buf,sizeof(buf),1,stdin); // read a word!
  while (success) {
    if ((buf & 0xffff) == GZ) { // if found GZ
      if((iWordInBlock-2 != nWordsInBlock) && nWordsInBlock)  std::cout<<"Warning::Found GZ in data"<<std::endl;
      else{ // If GZ is found in the header, this is a new block
        iWordInBlock = 0;
        nWordsInBlock = DecodeNDataInBlock(buf);
        if ((buf & 0x80ffff0000L) == START && len) { // if it's a new fragment, fill the tree.
          tree.Fill();
          len = 0;
	  vec.clear();
	  n++;
        }
      }
    }
    
    if (len == MAX_WORDS) {
      std::cout << "Too many words in fragment! (" << MAX_WORDS << " 64-bit words)" << std::endl;
      return 1;
    }
    
    if(SaveHeader || (iWordInBlock > 1)){//choose to store the block header
      blob[len] = buf; // store the read word
      vec.push_back(buf); // store the word into vector
      len++;
    }    
    iWordInBlock++;
    success = fread(&buf,sizeof(buf),1,stdin); // read another word!
  }
  tree.Fill();
  tree.Write();
  a.Close();
  return 0;
}
