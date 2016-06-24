(* Package written by 
      Pekka Janhunen
      Finnish Meteorological Institute
      Geophysics Dept.
      pjanhune@finsun.csc.fi  *)

BeginPackage["EinsteinTensor`"]

EinsteinTensor::usage = "EinsteinTensor[g,x] with g a nxn-matrix
  (the metric with lower indices) and x n-vector (the coordinates)
  gives the Einstein tensor (a nxn-matrix) with lower indices."

Begin["`Private`"]

EinsteinTensor[metric_,x_]:=
  Block[ {Dim,Metric, PreChristoffel, Christoffel, Riemann,
          PreRiemann, Ricci, CurvatureScalar,
          sigma, mu, nu, alpha, beta, gamma},
          Dim = Length[x];
          Metric = Simplify[Inverse[metric]];
            (* Metric with upper indices *)
          PreChristoffel =
            Table[ D[metric[[gamma,alpha]],x[[beta]]]
                 + D[metric[[beta,gamma]],x[[alpha]]]
           		    - D[metric[[alpha,beta]],x[[gamma]]],
           	   	 {gamma,Dim}, {alpha,Dim}, {beta,Dim} ];
           	 (* The "lower index part" of Christoffel symbols *)
          PreChristoffel = Simplify[PreChristoffel];
          Christoffel = (1/2) Metric . PreChristoffel;
             (* The full Christoffel symbols *)
          Christoffel = Simplify[Christoffel];
          PreRiemann = 
             Table[ D[Christoffel[[sigma,alpha,nu]],x[[mu]]]
                    + Sum[Christoffel[[gamma,alpha,nu]]
                            Christoffel[[sigma,gamma,mu]],
                          {gamma,Dim} ],
                    {sigma,Dim}, {alpha,Dim}, {mu,Dim}, {nu,Dim} ];
           	(* PreRiemann has to be antisymmetrized to yield
           	   Riemann tensor: *)
          Riemann = Table[ PreRiemann[[sigma,alpha,mu,nu]]
                         - PreRiemann[[sigma,alpha,nu,mu]],
                           {sigma,Dim}, {alpha,Dim},
                           {mu,Dim}, {nu,Dim} ];
          Ricci = Table[ Sum[Riemann[[sigma,alpha,sigma,beta]],
                             {sigma,Dim}],
                         {alpha,Dim}, {beta,Dim} ];
          CurvatureScalar = Sum[ Metric[[alpha,beta]]
                                 Ricci[[alpha,beta]],
                                 {alpha,Dim}, {beta,Dim} ];
          (* Return Einstein tensor: *)
          Ricci - (1/2) CurvatureScalar metric ]

End[]

EndPackage[]

Print[{EinsteinTensor}]

