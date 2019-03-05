// 
// Copyright(c) 1993-1996 Tony Kirke
// author="Tony Kirke" *
/* Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  This software is provided "as is" 
 * without express or implied warranty.
*/
#ifndef IIR1STHP
#define IIR1STHP
namespace SPUC {
//! \brief   Template for 1st Order IIR filter.
//
//!   The filter is assumed the first order low pass digital filter 
//!   generated by the bilinear transform of a first order analog 
//!   filter.
//!   G(z) =  (1+1/z)/(1-a/z) where a is real and |a| < 1
//!   Note: Not normalized
template <class Numeric> class iir_hp1
{
    protected:   
    	double gain;                    
    	double coeff;                    
    	Numeric out;
        Numeric previous_out;
		Numeric previous_in;
        
    public:
        iir_hp1(double A=0) {
		  set_coeff(A);
		  previous_in = previous_out = out = 0 ; 
		}
		void set_coeff(double A) { 
		  coeff=A;
		  gain=(1+coeff)/2.0;
		}
		//! Constructor reading coefficient from a file.
		iir_hp1(const char* file)
		  {
			FILE *iirf = fopen(file,"r"); 
			fscanf(iirf,"%lf",&coeff);
			fclose(iirf);
			previous_in = previous_out = out = 0;
			set_coeff(coeff);
		  }             
		//! Print out coefficients
		void print() {printf("IIR Coefficient gain = %lf\n",gain);}
		//! Input new sample and calculate output
		Numeric clock(Numeric input) {
		  // Shift previous outputs and calculate new output */
		  out = coeff*previous_out + (input-previous_in);
		  previous_out = out;
		  previous_in = input;
		  return(gain*out);
		}
		//! Reset
		void reset() {
		  previous_in = previous_out = out = 0;
		}
};      
} // namespace SPUC
#endif