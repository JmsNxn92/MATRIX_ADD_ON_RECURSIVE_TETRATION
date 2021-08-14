/* This program is the matrix add-on to my previous program Abel_L.gp. It runs without needing Abel_L.gp, but borrows functions from it.
This program also includes the file INIT.dat which stores a 100 x 100 matrix of Taylor coefficients which are enough for 100 digit precision.
If the INIT.DAT file is corrupted or lost; run the command,

beta_init(100);

This will take approximately an hour or two to compile. If it takes too long; you can redownload the .zip file this code was packaged in.

This code is entirely written by James David Nixon; and is an attachment to his construction of the beta tetration.
This code is coupled with the paper The Limits Of A Family; Of Asymptotic Solutions to the Tetration Equation. 
It can be found on the author's arxiv page.
*/


/*
set series precision/numerical precision to 100
A hundred isn't really needed, but for the sake of showing this converges for arbitrary precision, 100 is chosen.
To make it work for larger precision would require running and creating a larger INIT.dat.
This means you'd have to run beta_init(X) for X>100; I can't imagine how long that'd take though.
*/
\p 100
\ps 100

/*
A brief note on the function protocol. Everywhere you see a {v=0} is a marker for grabbing taylor series or not.
For example, if you were to write, Abel(1+z,log(2),z), this will print out the Taylor series about the point 1.
Or beta(3+z,log(2),z) will print out the Taylor series about the point 3. It's fairly self explanatory.
This code is not really based off of Taylor series; so it runs much slower than just numerical evaluation.
*/

/*******************************************************************************************/

/*
This is a construction of a faster running beta function, but it requires more data; and avoids much of the errors.
This program should have come with the file INIT.dat which includes all the Taylor series Coefficients needed for 100 digit accuracy.
To recreate this file, simply run beta_init(100); or beta_init(x) for any larger x; but 100 is enough for 100 digits of accuracy. 
I've coded everything for 100 digit accuracy or less; and if you want to meddle around there could be errors
Specifically, they'd be errors in the out = exp(out) process; where the exponential overflows.
We first use the function Phi_Inv to discover a taylorseries for beta(log(w)) about w=0, and then use this Taylor series and the functional equation to extend beta.
*/

Phi_Inv(w,l,{n=100}) =
{
	my(out = 0);
	for(i=0,n,
		out = w*exp(out)/(exp(l*(n+1-i))+w)
	);
	out;
}

/* this is the initialization protocol; you do not need to run it if you already have the INIT.dat file. */


beta_init(n) = {
	beta_taylor = Phi_Inv(w,l,n);
	writebin("INIT.DAT",beta_taylor);
	print("done.");
}

/*This just reads INIT.DAT; just make sure it's in the same folder as the source. And pari-gp is reading from the source's folder.*/

beta_taylor = read("INIT.DAT");

/*
This is the  beta function we need. 
It loads the Taylor values for small values, and then we push forward using the functional equation.
*/

beta(z,y,{v=0}) = {
	if(v==0,
		if(real(z) < -50, 
			sum(j=1,99,subst(Pol(polcoef(beta_taylor,j,w),l),l,y)*exp(y*z*j)),
			exp(beta(z-1,y,v))/(exp(-y*(z-1)) + 1)
		),
		if(real(polcoef(z,0,v)) < -50, 
			sum(j=1,99,subst(Pol(polcoef(beta_taylor,j,w),l),l,y)*exp(y*z*j)),
			exp(beta(z-1,y,v))/(exp(-y*(z-1)) + 1)
		)
	);
}

/*
This is the error function between beta and the tetration function associated with the multiplier y.
The value 3 is chosen because it is the largest value I could find which didn't cause fractal errors near the real line.
Luckily, it managed to be enough to get 100 digit precision.
Mathematically 3 can be any value; and in reality should be as large as possible, but in the confines of pari; we have to limit ourselves for overflows.
Ipso, making 3 larger usually results in fractal "hairs" on the produced function.
It'll look a little off and fractally; instead of being smooth and holomorphic.
Nonetheless, this still allows for 100 digit precision.
It limits us from achieving 1000 digit precision though; without producing fractal anomalies. 
Locally; everything will work unless we sit on a hair, though.
*/

rho(z,y,{v=0},{count=15})={
	if(v==0,
		if(abs(beta(z,y)) <= 3 && count>0, 
			count--;
			log(1 + (rho(z+1,y,v,count)-log(1+exp(-y*(z+1))))/beta(z+1,y)),
			0
		),
		if(abs(polcoef(beta(z,y,v),0,v)) <= 3 && count>0,
			count--;
			log(1 + (rho(z+1,y,v,count)-log(1+exp(-y*(z+1))))/beta(z+1,y)),
			0
		)
	);
}

/*
This is the first Abel function we construct.
It is associated to the multiplier y; which we will vary in the final result.

This Abel function will not be normalized as Abel(0,y) = 1; in order to do this, search for the zero at Abel(x_0,y) = 0 and redefine Abel(z,y) = Abel(z+2+x_0,y).
To do this, just convert the Taylor series into a polynomial and use pari's built in root finder.
*/

Abel(z,y,{v=0}) = {
	if(v==0,
		if(real(z) <= 0,
			beta(z,y,v) + rho(z,y,v) - log(1+exp(-y*z)),
			exp(Abel(z-1,y,v))
		),
		if(real(polcoef(z,0,v)) <=0,
			beta(z,y,v) + rho(z,y,v)-log(1+exp(-y*z)),
			exp(Abel(z-1,y,v))
		),
	);
}

/*This is the warped beta function which will produce the correct tetration*/

beta_warped(z,{v=0}) = {
	if(z==-1,z=0); /*catches a pointwise error*/
	if(v==0,
		beta(z,1/sqrt(z+1)),
		beta(z,1/sqrt(z+1),v)
	);
}

/*
This is the warped error function between beta_warped and the correct tetration.
This function is less than optimal, but will work; expect anomalies though.
The value 4 is chosen similarly to how 3 was chosen before. Choosing a larger constant will cause overflows very fast.
Especially with beta_warped; which is more anomalous than beta itself. It is volatile.
This is enough to get 100 digit precision though; excepting anomalous points.

To handle anomalous points; I've included two different methods.

***use the Taylor expansion protocol. 
***use the rho method, which is superior.
*/

tau_warped(z,{v=0},{count=8})={
	if(v==0,
		if(abs(beta_warped(z)) < 4 && count>0, 
			count--; 
			log(1 + tau_warped(z+1,v,count)/beta_warped(z+1)) +beta(z,1/sqrt(2+z)) - beta_warped(z) - log(1+exp(-z/sqrt(2+z))),
			0
		),
		if(abs(polcoef(beta_warped(z,v),0,v))<4 && count>0,
			count--;
			log(1 + tau_warped(z+1,v,count)/beta_warped(z+1,v)) +beta(z,1/sqrt(2+z),v) - beta_warped(z,v) - log(1+exp(-z/sqrt(2+z))),
			0
		)
	);
}

/*
This is the rho method; which takes a second order expansion (tau is first order).
This is difficult to explain without the user reading the paper paired with this code.
But, it ultimately results in faster running code; it just nests a composition rather than not doing so--like with tau_warped.
It's little more than a change of variables.  With infinite memory and infinite time, \tau_warped is perfectly coded--but life is life.

count can be increased to as much as you want; similarly with the limiter 10; but upon which, expect a plethora of overflow errors.
Or expect fractal anomalies.
I've set count to 30 by default; as this gives near perfect results.

There are still anomalies with the rho method. Grabbing Taylor series and summing them is still the best method.
*/

rho_warped(z,{v=0},{count=30})={
	if(v==0,
		if(abs(beta_warped(z)) < 3 && count>0, 
			count--; 
			log((rho_warped(z+1,v,count) +beta(z+1,1/sqrt(3+z)) - log(1+exp(-(z+1)/sqrt(3+z))))/beta_warped(z+1)),
			0
		),
		if(abs(polcoef(beta_warped(z,v),0,v)) < 3 && count>0,
			count--;
			log((rho_warped(z+1,v,count) +beta(z+1,1/sqrt(3+z),v) - log(1+exp(-(z+1)/sqrt(3+z))))/beta_warped(z+1,v)),
			0
		)
	);
}


/*
This is the correct super exponential. First the not normalized one, and the normalized one. 
The normalization constant x_0 = 1.969637739698065306544624079350257708852542229771084623924562193889980567396859073585886113019833431
is calculated by x_0 = 2 + polrootsreal(Pol(rho_warped(1+z,z) +beta(z+1,1/sqrt(3+z),z) - log(1+exp((-z-1)/sqrt(3+z))),z))
If you fiddle with the code, the first thing that tends to stray is the normalization constant, so just run that command and you'll get it back.
*/

x_0 = 1.969637739698065306544624079350257708852542229771084623924562193889980567396859073585886113019833431;

Sexp_NOT_NORMALIZED(z,{v=0},{count=30}) = {
	if(v==0,
		if(real(z) <= 1,
			rho_warped(z,v,count) +beta(z,1/sqrt(2+z)) - log(1+exp(-z/sqrt(2+z))),
			exp(Sexp_ NOT_NORMALIZED(z-1,v))
		),
		if(real(polcoef(z,0,v)) <= 1,
			rho_warped(z,v,count) +beta(z,1/sqrt(2+z),v) - log(1+exp(-z/sqrt(2+z))),
			exp(Sexp_NOT_NORMALIZED(z-1,v))
		)
	);
}

/*************************************************
DISCLAIMER:

This is the final tetration function. It still produces errors and anomalous points because it's based off the rho method.
To get a more accurate picture; the author suggests expanding a taylor series near the anomalous points.
This would be Sexp(A + z,z) to get the taylor series at A; the command Sexp(A + z) will not run. For A a number.
Be sure to flag the variable you are taking the taylor series in.
Remember the radius of convergence.

Please remember, this is only rough code.
It's beta code of the beta method ;D
**************************************************/

Sexp(z,{v=0}, {count=30}) = {
	Sexp_NOT_NORMALIZED(z+x_0,v,count);
}


/* HERE ENDS ALL THE CODE I'VE WRITTEN--the rest was procured*/


/* =============================================================================  
** Color graphing system - mike3 - 20.11.10 04:07
** Hi.
** I thought I'd post the code I use to generate the color graphs from Pari/GP.
** Here it is.
**
** Note: the output is in .PPM format. 
** You'll need something else to convert that to .PNG. (I use GIMP.)
** 
** Also, I might warn you: it takes a LONG time to graph a complex function
** with a significantly complicated calculation procedure, as it must be
** evaluated at every pixel of the graph.
** 
** (updated 12/16/2010 -- program was not all good: 
**   * spurious "func" parameter in MakeGraph and "safetyarg" was missing.)
** ------------------------------------------------------------------------------------------ 
**  
** ============================================================================= */

/* =============================== Code: ==================================================== */
/* Complex function magnitude/phase plotter. */

/* To use:
*     1. Define function to graph as func(z).
*     2. Load this program.
*     3. Execute MakeGraph(width, height, x0, y0, x1, y1, filename) with the parameters given as follows:
*        width, height = width/height of image in pixels
*        x0, y0, x1, y1 = rectangle of complex plane to graph: x0 + y0i in upper-left corner to x1 + y1i in lower-right corner
*        filename = name of file to save as.
* Output is in .PPM format.
*/


/* Color conversion (HSB to RGB). */

HSB2RGB(hsb) = {
       local(H=hsb[1]);
       local(S=hsb[2]);
       local(B=hsb[3]);
       local(HH);
       local(F);
       local(P);
       local(Q);
       local(T);

       HH = floor(6*H)%6;
       F = (6*H) - floor(6*H);
       P = B*(1 - S);
       Q = B*(1 - (F*S));
       T = B*(1 - (1-F)*S);
       if(B > 1.0, B = 1.0);
       if(B < 0.0, B = 0.0);
       if(P > 1.0, P = 1.0);
       if(P < 0.0, P = 0.0);
       if(Q > 1.0, Q = 1.0);
       if(Q < 0.0, Q = 0.0);
       if(T > 1.0, T = 1.0);
       if(T < 0.0, T = 0.0);

       if(HH == 0, return([B, T, P]));
       if(HH == 1, return([Q, B, P]));
       if(HH == 2, return([P, B, T]));
       if(HH == 3, return([P, Q, B]));
       if(HH == 4, return([T, P, B]));
       if(HH == 5, return([B, P, Q]));
       }

/* Safe argument. */
safetyarg(z) = if(z == 0, 0, arg(z));

/* Make graph. */
MakeGraph(width, height, x0, y0, x1, y1, filename) = {
       xstep = (x1 - x0)/width;
       ystep = (y1 - y0)/height;
       write(filename, "P3");
       write(filename, "# ", filename);
       write(filename, width, " ", height);
       write(filename, "255");

       for(y=0, height-1,
           for(x=0, width-1,
                  xx = x0+(xstep*x);
                  yy = y0+(ystep*y);
               z = xx+yy*I;
               funcvalue = func(z);
               mag = abs(funcvalue);
               phase = safetyarg(funcvalue);
               H = phase/(2*Pi);
               S = 1/(1 + 0.3*log(mag + 1));
               B = 1 - 1/(1.1 + 5*log(mag + 1));
               RGB = HSB2RGB([H, S, B]);
                  Red = floor(RGB[1]*255.0);
                  Green = floor(RGB[2]*255.0);
                  Blue = floor(RGB[3]*255.0);
               write1(filename, Red, " ", Green, " ", Blue, "  ");
              );
           write(filename, "");
       );
       print("Done.");
    } 
