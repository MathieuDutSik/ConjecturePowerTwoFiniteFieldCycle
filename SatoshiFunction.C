#include <NTL/ZZX.h>
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>

#define nDim 25

NTL_CLIENT
// one function for mapping 01 vector to integer
// one function for given a 01 vector, find the next one.
// one function that takes the integer h and return the list of x^(2^h-1)
// one list of mapping x --> x+1
// one function that takes two prime integers h, k and returns
//   the cycle decomposition

void pow2Func(long & h, long & pow)
{
  int i;
  pow=1;
  for (i=1; i<=h; i++)
    {
      pow=pow*2;
    }
}


void PrintSequence(vec_long & OneSeq)
{
  long i;
  cout << "(";
  for (i=0; i<OneSeq.length(); i++)
    {
      if (i>0)
	{
	  cout << ",";
	}
      cout << OneSeq(i);
    }
  cout << ");\n";
}


int GCD(long & a, long & b)
{
  long aP, bP;
  aP=a;
  bP=b;
  while( 1 )
    {
      aP = aP % bP;
      if( aP == 0 )
	return bP;
      bP = bP % aP;
      if( bP == 0 )
	return aP;
    }
}


// returns x, y such that xa+yb=gcd(a,b)
void BezoutTheorem(long & a, long & b, long & x, long & y)
{
  long q, r, xN, yN;
  if (b == 0)
    {
      x=1;
      y=0;
    }
  else
    {
      r=a%b;
      q=(a-r)/b;
      BezoutTheorem(b,r,xN, yN);
      x=yN;
      y=xN-q*yN;
    }
}



void Mapping01toLong(GF2X& The01, long&  TheZint)
{
  long i;
  long pow;
  long OtherPow;
  long prov;
  GF2 Pre;
  prov=nDim-1;
  TheZint=0;
  OtherPow=1;
  for (i=0; i<=prov; i++)
    {
      Pre=coeff(The01, i);
      TheZint=TheZint+rep(Pre)*OtherPow;
      OtherPow=OtherPow*2;
    }
}







void NextFieldElement(GF2X& TheElt)
{
  long place;
  long testFinish;
  long i;
  place=0;
  while (rep(coeff(TheElt, place))==1)
  {
    place++;
    if (place>=nDim)
      {
	testFinish=1;
	break;
      }
  }
  if (testFinish==1)
    {
      for (i=0; i<nDim; i++)
	{
	  SetCoeff(TheElt, i, 0);
	}
    }
  else
    {
      for (i=0; i<place; i++)
	{
	  SetCoeff(TheElt, i, 0);
	}
      SetCoeff(TheElt, i, 1);
    }
}





// this function returns a generating element 
// of the multiplicative group of GF(2^nDim)
void GetAGeneratingElement(GF2X & x, GF2XModulus & IrrP)
{
  int IsOK;
  long prov, i, pow, pos, posDivisor;
  GF2X NeutralElt, ThePow;
  long re;
  vec_long ListDivisor;
  x=0;
  NeutralElt=0;
  SetCoeff(NeutralElt, 0, 1);
  prov=nDim;
  pow2Func(prov, pow);
  pow=pow-1;
  ListDivisor.SetLength(pow);
  posDivisor=1;
  for (i=1; i<pow; i++)
    {
      re=pow%i;
      if (re == 0)
	{
	  ListDivisor(posDivisor)=i;
	  posDivisor++;
	}
    }
  while(1)
    {
      NextFieldElement(x);
      IsOK=1;
      pos=1;
      while(1)
	{
	  ThePow=PowerMod(x, ListDivisor(pos), IrrP);
	  if (ThePow == NeutralElt)
	    {
	      IsOK=0;
	      break;
	    }
	  pos++;
	  if (pos == posDivisor)
	    {
	      break;
	    }
	}
      if (IsOK == 1)
	{
	  break;
	}
    }
}


void PowerSequence(GF2XModulus & IrrP, vec_long & PowerSeq, vec_long & RevPowerSeq, GF2X & TheGenerator)
{
  long prov, pow, i, pos;
  GF2X SuccPow, provGF2X;
  long ePos;
  prov=nDim;
  pow2Func(prov, pow);
  pow=pow-1;
  PowerSeq.SetLength(pow);
  RevPowerSeq.SetLength(pow);
  SetCoeff(SuccPow, 0, 1);
  for (i=0; i<pow; i++)
    {
      Mapping01toLong(SuccPow, ePos);
      PowerSeq(ePos-1)=i;
      RevPowerSeq(i)=ePos-1;
      mul(provGF2X, SuccPow, TheGenerator);
      rem(SuccPow, provGF2X, IrrP);
    }
  //  PrintSequence(PowerSeq);
  //  PrintSequence(RevPowerSeq);
}



void TheSequenceVers2(long & h, long & k, vec_long & PowerSeq, vec_long & RevPowerSeq, vec_long & RetSeq)
{
  long prov, pow, pos, i;
  long powD, powH, powK, x, y, fact;
  ZZ provZZ, powDZZ;
  prov=nDim;
  pow2Func(prov, pow);
  powD=pow-1;

  pow2Func(h, pow);
  powH=pow-1;
  
  pow2Func(k, pow);
  powK=pow-1;
  
  BezoutTheorem(powD, powH, x, y);
  ZZ prov1, prov2, prov3;
  prov1=to_ZZ(y);
  prov2=to_ZZ(powK);
  mul(prov3, prov1, prov2);
  prov1=to_ZZ(powD);
  QuickRem(prov3, prov1);
  fact=to_long(prov3);

  pow2Func(prov, pow);
  RetSeq.SetLength(pow);
  powDZZ=to_ZZ(powD);
  for (i=0; i<pow; i++)
    {
      if (i == 0)
	{
	  RetSeq(i)=0;
	}
      else
	{
	  provZZ=to_ZZ(PowerSeq(i-1))*to_ZZ(fact);
	  QuickRem(provZZ, powDZZ);
	  RetSeq(i)=RevPowerSeq(to_long(provZZ))+1;
	}
    }


}




void SequenceOfAdd1(vec_long & TheSeq)
{
  long i, pow, prov;
  long eVal;
  prov=nDim;
  pow2Func(prov, pow);
  TheSeq.SetLength(pow);
  prov=nDim-1;
  pow2Func(prov, pow);
  for (i=0; i<pow; i++)
    {
      TheSeq(2*i)=2*i+1;
      TheSeq(2*i+1)=2*i;
    }
}


void ReverseSeq(vec_long & TheSeq, vec_long & RevSeq)
{
  long pos, pos2, i, pow;
  long posZZ;
  pow=1;
  for (i=0; i<nDim; i++)
    {
      pow=pow*2;
    }
  RevSeq.SetLength(pow);
  pow=pow-1;
  pos=0;
  while(1)
    {
      posZZ=pos;
      pos2=0;
      while(1)
	{
	  if (TheSeq(pos2) == posZZ)
	    {
	      break;
	    }
	  pos2++;
	}
      RevSeq(pos)=pos2;
      if (pos == pow)
	{
	  break;
	}
      pos++;
    }
}



void ComposeSequence(vec_long & Seq1, vec_long & Seq2, vec_long & NewSeq)
{
  long pow, i, prov;
  long eVal;
  prov=nDim;
  pow2Func(prov, pow);
  for (i=0; i<pow; i++)
    {
      eVal=Seq2(i);
      NewSeq(i)=Seq1(eVal);
    }
}






void SatoshiSequenceVers2(long& h, long &k, vec_long & PowerSeq, vec_long & RevPowerSeq, vec_long & TheReply)
{
  vec_long Add1;
  vec_long Epsilon, InvEpsilon;
  vec_long Prov;
  long prov, pow;
  prov=nDim;
  pow2Func(prov, pow);
  Epsilon.SetLength(pow);
  InvEpsilon.SetLength(pow);
  Prov.SetLength(pow);
  TheSequenceVers2(h, k, PowerSeq, RevPowerSeq, Epsilon);
  TheSequenceVers2(k, h, PowerSeq, RevPowerSeq, InvEpsilon);
  SequenceOfAdd1(Add1);
  ComposeSequence(Add1, Epsilon, TheReply);
  ComposeSequence(InvEpsilon, TheReply, Prov);
  ComposeSequence(Add1, Prov, TheReply);

}





void ExtractTheCycle(vec_long & TheSeq, vec_long & CycleSequence)
{
  vec_long Status;
  vec_long ListLengthCycle;
  vec_long Magnitude;
  long pos, WorkPos, prov, pow;
  long nbCycle, len, i, Beginpos;
  int IsFinished, eVal, eMag;
  long eTest;
  eTest=1;
  prov=nDim;
  pow2Func(prov, pow);
  Status.SetLength(pow);
  Magnitude.SetLength(pow);
  CycleSequence.SetLength(pow);
  for (i=0; i<pow; i++)
    {
      Status(i)=1;
      CycleSequence(i)=0;
      Magnitude(i)=0;
    }
  nbCycle=0;
  while(1)
    {
      pos=0;
      IsFinished=0;
      while(1)
	{
	  if (Status(pos)==eTest)
	    {
	      Beginpos=pos;
	      break;
	    }
	  pos++;
	  if (pos==pow)
	    {
	      IsFinished=1;
	      break;
	    }
	}
      if (IsFinished==1)
	{
	  break;
	}
      WorkPos=Beginpos;
      len=0;
      while(1)
	{
	  Status(WorkPos)=0;
	  len++;
	  WorkPos=TheSeq(WorkPos);
	  if (WorkPos == pos)
	    {
	      break;
	    }
	}
      CycleSequence(nbCycle)=len;
      nbCycle++;
    }
  for (i=0; i<pow; i++)
    {
      eVal=CycleSequence(i);
      if (eVal >0)
	{
	  Magnitude(eVal)=Magnitude(eVal)+1;
	}
    }
  for (i=0; i<pow; i++)
    {
      eMag=Magnitude(i);
      if (eMag>0)
	{
	  cout << "  " << i << "^" << eMag;
	}
    }
}





void ExtractTheCycleVers2(vec_long & TheSeq, vec_long & CycleSequence)
{
  vec_long ListLengthCycle;
  vec_long Magnitude;
  vec_long ListNext, ListPrev;
  long pos, WorkPos, BeginPos, prov, pow;
  long nbCycle, len, i;
  int eVal, eMag, eNext, ePrev;
  long eTest;
  eTest=1;
  prov=nDim;
  pow2Func(prov, pow);
  Magnitude.SetLength(pow);
  CycleSequence.SetLength(pow);
  ListNext.SetLength(pow+2);
  ListPrev.SetLength(pow+2);
  ListNext(pow)=0;
  ListPrev(pow+1)=pow-1;
  for (i=0; i<pow; i++)
    {
      CycleSequence(i)=0;
      Magnitude(i)=0;
      if (i<pow-1)
	{
	  ListNext(i)=i+1;
	}
      else
	{
	  ListNext(i)=pow+1;
	}
      if (i == 0)
	{
	  ListPrev(i)=pow;
	}
      else
	{
	  ListPrev(i)=i-1;
	}
    }
  nbCycle=0;
  while(1)
    {
      eNext=ListNext(pow);
      if (eNext == pow+1)
	{
	  break;
	}
      BeginPos=eNext;
      WorkPos=eNext;
      len=0;
      while(1)
	{
	  ePrev=ListPrev(WorkPos);
	  eNext=ListNext(WorkPos);
	  ListNext(ePrev)=eNext;
	  ListPrev(eNext)=ePrev;
	  len++;
	  WorkPos=TheSeq(WorkPos);
	  if (WorkPos == BeginPos)
	    {
	      break;
	    }
	}
      //      cout << "Before\n";
      CycleSequence(nbCycle)=len;
      //      cout << "After\n";
      nbCycle++;
    }
  for (i=0; i<pow; i++)
    {
      eVal=CycleSequence(i);
      if (eVal >0)
	{
	  Magnitude(eVal)=Magnitude(eVal)+1;
	}
    }
  for (i=0; i<pow; i++)
    {
      eMag=Magnitude(i);
      if (eMag>0)
	{
	  cout << "  " << i << "^" << eMag;
	}
    }
}



int main()
{
  GF2X P, TheGen;
  GF2XModulus Pmod;
  vec_long TheSequence, PowerSeq, RevPowerSeq;
  vec_long CycleSeq;
  long i;
  long h, k, nM;
  long prov, pow;
  BuildSparseIrred(P, nDim);
  build(Pmod, P);
  nM=nDim-1;
  GetAGeneratingElement(TheGen, Pmod);
  PowerSequence(Pmod, PowerSeq, RevPowerSeq, TheGen);
  prov=nDim;
  pow2Func(prov, pow);
  TheSequence.SetLength(pow);
  for (h=1; h<nDim; h++)
    {
      if (GCD(h,nDim) == 1)
	{
	  for (k=h+1; k<nDim; k++)
	    {
	      if (GCD(k,nDim) == 1)
		{
		  SatoshiSequenceVers2(h, k, PowerSeq, RevPowerSeq, TheSequence);
		  cout << "h=" << h << " k=" << k << ":  ";
		  // ExtractTheCycle(TheSequence, CycleSeq);
		  ExtractTheCycleVers2(TheSequence, CycleSeq);
		  cout << "\n";
		}
	    }
	}
    }
}
