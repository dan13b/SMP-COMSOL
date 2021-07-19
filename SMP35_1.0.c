/* Interface to a Shape Memory Polymer 	*/

/* Version 1.0 created on 19-07-2021  */
/* Daniel Bautista-Salinas       		*/


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

typedef struct EXPORT{
 float m[6][6];
} arr2d;

EXPORT double determinant(double matrix[6][6], double size);

EXPORT arr2d cofactor(double matrix[6][6], double size);

EXPORT arr2d transpose(double matrix[6][6], double matrix_cofactor[6][6], double size);


EXPORT int eval(double e[6],        // Input: Green-Lagrange strain tensor components in Voigt order (xx,yy,zz,yz,zx,xy)
                double s[6],        // Output: Second Piola-Kirchhoff stress components in Voigt order (xx,yy,zz,yz,zx,xy)
                double Cei[6][6],   // Output: Jacobian of stress with respect to strain, 6-by-6 matrix in row-major order
                int *nPar,        	// Input: Number of material model parameters, scalar
                double *par,        // Input: Parameters: par[0] = E, par[1] = nu
                int *nStates1,      // Input: Number of states1 (scalar)
                double *states1,    // States: Volume fraction of the glassy phase, xi_g
				int *nStates2,		// Input: Number of states2 (scalar)
				double *states2) {  // States: Inelastic strain, eps_i

	int i,j,k;
	int size=6;
	double d;

    double Sr[6][6];
    double Sg[6][6];
    double c_e[6][6];
    int I[6][6];

	double alpha_r, alpha_g, E_r, E_g, nu_r, nu_g;
	double a, b, T, Tl,Th, Tg, Tr;
	double x;

	double xi_g;
	double deltaXi;
	double eps_t;
    double eps_i[6];
    double eps[6];
    double mu;

    double A[6][6];
	double B[6][6];
	double C[6][6];
	double D[6][6];
	double E[6][6];
	double F[6];
	double sum = 0.0;

	arr2d Ce;
	arr2d Ainv;
	arr2d Cinv;
	arr2d Dinv;

	// Check inputs
	if (nPar[0] != 14)           // Number of parameters needed
	return 1;                    // error code 1 = "Wrong number of parameters"
	if (nStates1[0] !=1)         // Number of states needed
	return 2;                    // error code 2 = "Wrong number of states"
	if (nStates2[0] !=6)         // Number of states needed
	return 2;					 // error code 2 = "Wrong number of states"

	// Read input parameters from COMSOL
	alpha_r = par[0];
	alpha_g = par[1];
	E_r = par[2];
	E_g = par[3];
	nu_r = par[4];
	nu_g = par[5];
	a = par[6];
	b = par[7];
	T = par[8];
	Tl = par[9];
	Th = par[10];
	Tg = par[11];
	Tr = par[12];
	x = par[13];

	// Check values of input parameters (To be completed)
	if (alpha_r != 0.0001) return -1;
	if (alpha_g != 0.0001) return -1;
	if (E_r <= 0.0) return -1;
	if (E_g <= 0.0) return -1;
	if (nu_r >= 0.5 || nu_r <= -1.0) return -1;
	if (nu_g >= 0.5 || nu_g <= -1.0) return -1;
	if (a != 0.185) return -1;
	if (b != 0.180) return -1;
	if (T > Th || T < Tl) return -1;
	if (Tl != 278.15) return -1;
	if (Th != 358.15) return -1;
	if (Tg != 308.15) return -1;
	if (Tr != 297.15) return -1;
	if (x > 5.0) return -1;

	
	// 6x6 compliance matrix of glassy and rubbery phases and overall equivalent stiffness. Initially filled with zeros  
	for (i = 0; i < 6; i++){
		for (j = 0; j < 6; j++){
	    	Sr[i][j] = 0.0;
	    	Sg[i][j] = 0.0;
            c_e[i][j] = 0.0;
            Cei[i][j] = 0.0;
            I[i][j] = 0.0; // identity matrix
            A[i][j] = 0.0;
			B[i][j] = 0.0;
			C[i][j] = 0.0;
			D[i][j] = 0.0;
			E[i][j] = 0.0;
		}
		eps_i[i] = 0;
	}

	// Read state varaibles and initial calculations
    for (i=0;i<6;i++)  eps_i[i] = states2[i];

	Sr[0][0] = Sr[1][1] = Sr[2][2] = 1/E_r; // upper diagonal
	Sr[3][3] = Sr[4][4] = Sr[5][5] = 2*(1+nu_r)/E_r; // lower diagonal
	Sr[0][1] = Sr[0][2] = Sr[1][0] = Sr[1][2] = Sr[2][0] = Sr[2][1] = -nu_r/E_r; // off diagonal, equals lambda Lame

	Sg[0][0] = Sg[1][1] = Sg[2][2] = 1/E_g; // upper diagonal
	Sg[3][3] = Sg[4][4] = Sg[5][5] = 2*(1+nu_g)/E_g; // lower diagonal
	Sg[0][1] = Sg[0][2] = Sg[1][0] = Sg[1][2] = Sg[2][0] = Sg[2][1] = -nu_g/E_g; // off diagonal, equals lambda Lame

	I[0][0] = I[1][1] = I[2][2] = I[3][3] = I[4][4] = I[5][5] = 1;

	xi_g = -(tanh(a*(Tg-273.15)-b*(T-273.15))-tanh(a*(Tg-273.15)-b*(Th-273.15)))/(tanh(a*(Tg-273.15)-b*(Th-273.15))-tanh(a*(Tg-273.15)-b*(Tl-273.15))); // Volume fraction of the glassy phase

	eps_t = alpha_r*(T-Tr); // Thermal strain

 	for (i = 0; i < 6; i++){
		for (j = 0; j < 6; j++){
   			c_e[i][j] = (Sr[i][j]+xi_g*(Sg[i][j]-Sr[i][j])); // Overall equivalent stiffness
		}
	}

	d=determinant(c_e,size);
	if (d==0){
		return -1;
	}
	else
		Ce = cofactor(c_e,size);

	deltaXi = xi_g - states1[0];


	// x is the parameter used to move between the different SMP states
	// x<1 is the heating stage between T_low (x=0) and T_high (x=1)
	// x >= 1 && x <= 2.0 is the loading stage at a constant temperature, T_high
	// x > 2.0 && x < 3.0 is the cooling stage between T_high (x=2) and T_low (x=3)
	// x >= 3 && x <= 4 is the reheating/recovery stage between T_low (x=3) and T_high (x=4)


	if (x < 1.0){
		for(i = 0; i < 6; i++){
            for(j = 0; j < 6; j++){
                Cei[i][j] = Ce.m[i][j]; // Material Jacobian
            }
        }

	 	mu = 0;	

    	for (i = 0; i < 6; i++){
			eps[i] = eps_t; // Heating (1->2)
   	 	}

   	 	// Compute stress tensor: s = C*(e-et-ei)
		for (i = 0; i < 6; i++){
			s[i] = 0.0;
		}

		// Update states
		states1[0] = xi_g;
    	return 0; // Return value 0 if success, any other value triggers an exception
	}

	else if (x >= 1 && x <= 2.0){
        for(i = 0; i < 6; i++){
            for(j = 0; j < 6; j++){
                Cei[i][j] = Ce.m[i][j];
            }
        }

	 	mu = 0;	

    	for (i = 0; i < 6; i++){
			eps[i] = e[i]; // Loading (2->3)
   	 	}

   	 	// Compute stress tensor: s = C*(e-et-ei)
		for (i = 0; i < 6; i++){
			s[i] = 0.0;
			for (j = 0; j < 6; j++){
				s[i] += Cei[i][j]*(e[j]-eps_t-(mu*eps_i[j]));
			}
		}
		//for (i=0;i<6;i++)  eps_i[i] = mu*states2[i]; // Calculate new inelastic strain (it will be 0 during heating)

		// Update states
		states1[0] = xi_g;
        //for (i=0;i<6;i++)  states2[i] = eps_i[i];
    	return 0;
	}	
	
	else if (x > 2.0 && x < 3.0){
	 	// Cei during cooling
	 		// Calculate A
			for (i = 0; i < 6; i++){
			  	for (j = 0; j < 6; j++){
			    	for (k = 0; k < 6; k++){
			      		sum = sum + Sr[i][k]*Ce.m[k][j];
			    	}
			    A[i][j] = sum;
			    sum = 0;
			  	}
			}

			// Calculate Ainv
			d=determinant(A,size);
			if (d==0){
			    return -1;
			}
			else
			    Ainv = cofactor(A,size);

			// Calculate B
			for (i = 0; i < 6; i++){
			  	for (j = 0; j < 6; j++){
			  		B[i][j] = deltaXi*A[i][j];
				}
			}

			// Calculate C
			for (i = 0; i < 6; i++){
			  	for (j = 0; j < 6; j++){
			  		C[i][j] = I[i][j]+B[i][j];
				}
			}

			// Calculate Cinv
			d=determinant(C,size);
			if (d==0){
			    return -1;
			}
			else
			    Cinv = cofactor(C,size);

			// Calculate D
			for (i = 0; i < 6; i++){
			  	for (j = 0; j < 6; j++){
			  		D[i][j] = deltaXi*I[i][j] + Ainv.m[i][j];
				}
			}

			// Calculate Dinv
			d=determinant(D,size);
			if (d==0){
			    return -1;
			}
			else
			    Dinv = cofactor(D,size);

			// Calculate E
			for (i = 0; i < 6; i++){
			  	for (j = 0; j < 6; j++){
			  		E[i][j] = I[i][j] - deltaXi*Dinv.m[i][j];
				}
			}

			// Calculate Cei
			for (i = 0; i < 6; i++){
			  	for (j = 0; j < 6; j++){
			    	for (k = 0; k < 6; k++){
			      		sum = sum + Ce.m[i][k]*E[k][j];
			    	}
			    Cei[i][j] = sum;
			    sum = 0;
			  	}
			}

   	 	mu = 1;

	 	for (i = 0; i < 6; i++){
	 		// if (T>Tl){
	 			eps[i] = e[i]; // Cooling (3->4)
	 		// }
	 		// else{
	 		// 	eps[i] = eps_t + eps_i[i]; // Unloading (4->5)
	 		// }
	 	}

        for (i = 0; i < 6; i++){
            s[i] = 0.0;
            // if (T>Tl){
            	for (j = 0; j < 6; j++){
                	s[i] += Cei[i][j]*(e[j]-eps_t-eps_i[j]); //Stress (it will be 0 during unloading)
            	}
        	// }
        }

        // if (T>Tl){
			// Calculate F
			for (i = 0; i < 6; i++){
				F[i] = 0.0;
			  	for (j = 0; j < 6; j++){
			  		F[i] += B[i][j]*(eps[j]-eps_t);
				}
			  	F[i] += states2[i];
			}

			// Update eps_i
			for (i = 0; i < 6; i++){
				eps_i[i] = 0;
			  	for (j = 0; j < 6; j++){
			  		eps_i[i] += Cinv.m[i][j]*F[j];
				}
			}
		// }

		// Update states
		states1[0] = xi_g;
    	for (i=0;i<6;i++)  states2[i] = 0;
		for (i=0;i<6;i++)  states2[i] = eps_i[i];
    	return 0;		
	}

	else{
        for(i = 0; i < 6; i++){
            for(j = 0; j < 6; j++){
                Cei[i][j] = Ce.m[i][j];
            }
        }

	 	mu = xi_g/states1[0];

    	for (i = 0; i < 6; i++){
			eps[i] = eps_t + eps_i[i]; // Reheating (5->2)
   	 	}

   	 	// Compute stress tensor: s = C*(e-et-ei)
		for (i = 0; i < 6; i++){
			s[i] = 0.0;
			for (j = 0; j < 6; j++){
				s[i] += Cei[i][j]*(e[j]-eps_t-(mu*eps_i[j]));
			}
		}
		for (i=0;i<6;i++)  eps_i[i] = 0;
		for (i=0;i<6;i++)  eps_i[i] = mu*states2[i];

		// Update states
		states1[0] = xi_g;
		for (i=0;i<6;i++)  states2[i] = 0;
        for (i=0;i<6;i++)  states2[i] = eps_i[i];

    	return 0;
	}
}


// Matrix determinant (recursive function)
EXPORT double determinant(double matrix[6][6], double size){
    double s=1,det=0,m_minor[6][6];
    int i,j,m,n,c;
    if (size==1)
    {
        return (matrix[0][0]);
    }
    else
    {
        det=0;
        for (c=0;c<size;c++)
        {
            m=0;
            n=0;
            for (i=0;i<size;i++)
            {
                for (j=0;j<size;j++)
                {
                    m_minor[i][j]=0;
                    if (i != 0 && j != c)
                    {
                       m_minor[m][n]=matrix[i][j];
                       if (n<(size-2))
                          n++;
                       else
                       {
                           n=0;
                           m++;
                       }
                    }
                }
            }
            det=det + s * (matrix[0][c] * determinant(m_minor,size-1));
            s=-1 * s;
        }
    }
 
    return (det);
}
 
// Matrix cofactor
EXPORT arr2d cofactor(double matrix[6][6], double size){
    double m_cofactor[6][6],matrix_cofactor[6][6];
    arr2d invCe;
    int p,q,m,n,i,j;
    for (q=0;q<size;q++)
    {
        for (p=0;p<size;p++)
        {
            m=0;
            n=0;
            for (i=0;i<size;i++)
            {
                for (j=0;j<size;j++)
                {
                    if (i != q && j != p)
                    {
                        m_cofactor[m][n]=matrix[i][j];
                        if (n<(size-2))
                            n++;     
                        else
                        {
                            n=0;
                            m++;
                        }
                    }
                }
            }
        matrix_cofactor[q][p]=pow(-1,q + p) * determinant(m_cofactor,size-1);
        }
    }

    invCe = transpose(matrix,matrix_cofactor,size);
    return (invCe);
}

// Transpose of matrix cofactor
EXPORT arr2d transpose(double matrix[6][6], double matrix_cofactor[6][6], double size){
    int i,j;
    double m_transpose[6][6],d;
    arr2d m_inverse;

    for (i=0;i<size;i++)
    {
        for (j=0;j<size;j++)
        {
            m_transpose[i][j]=matrix_cofactor[j][i];
        }
    }
    d=determinant(matrix,size);
    for (i=0;i<size;i++)
    {
        for (j=0;j<size;j++)
        {
            m_inverse.m[i][j]=m_transpose[i][j] / d;
        }
    }
    return (m_inverse);
}