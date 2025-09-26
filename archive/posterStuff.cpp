


force = c/6*Ra*Ra*Ra*Rb*Rb*Rb*
	(h+Ra+Rb)/(hsq + twoRah + twoRbh)*
	(hsq + twoRah + twoRbh)*
	(hsq + twoRah + twoRbh + 4*Ra*Rb)*
	(hsq + twoRah + twoRbh + 4*Ra*Rb)


numerator = c*Ra*Ra*Ra*Rb*Rb*Rb*(h+Ra+Rb);
d1 = hsq + twoRah+ twoRbh;
d2 = d1 + 4*Ra*Rb;
denomRecip = 1/(6*d1*d1*d2*d2);
force = (numerator*denomRecip);



----------------------------------------------------

for(int A = 1; A < num_particles; +A)
{

	for (int B = 0; B < A; ++B)
	{

		// calculate forces
	}
}

----------------------------------------------------

len = num_particles;
#pragma omp parallel for
for (pc = 1; pc <= (((len*len)-len)/2); pc++)
{
    pd = (sqrt(pc*8.0+1.0)+1.0)*0.5;
    A = pd;
    B = (pd-A*(A-1.0)*.5-1.0);

    // calculate forces
}

----------------------------------------------------


\textrm{Iteration space } = \begin{bmatrix}
0 & {1,2} & {1,3} & {1,4} & {1,5} \\
0 & 0 & {2,3} & {2,4} & {2,5} \\
0 & 0 & 0 & {3,4} & {3,5} & \dots\\
0 & 0 & 0 & 0 & {4,5} \\
0 & 0 & 0 & 0 & 0 \\
&&\vdots&&&\ddots
\end{bmatrix}  