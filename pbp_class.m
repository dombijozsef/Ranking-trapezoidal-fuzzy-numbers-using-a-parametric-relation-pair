classdef pbp_class < handle
    
       
    methods(Static)    
        
        % This computes the preference index for two intervals 
        function M = M_interval(I_1,I_2)
            a1 = I_1(1,1);
            b1 = I_1(1,2);
            a2 = I_2(1,1);
            b2 = I_2(1,2);
            l1 = (b1-a1)/2;
            l2 = (b2-a2)/2;
            if l1>0 && l2>0
                d = (a2+b2)/2 - (a1+b1)/2;
                if a2<b2 && b2<=a1 && a1<b1
                    M = 0; 
                end
                if a1<b1 && b1<=a2 && a2<b2
                    M = 1; 
                end
                if a1<a2 && a2<b1 && b1<b2
                    M = 1-(l1+l2-d)^2/(8*l1*l2); 
                end
                if a2<a1 && a1<b2 && b2<b1
                    M = (l1+l2+d)^2/(8*l1*l2); 
                end
                if a2<=a1 && a1<b1 && b1<=b2
                    M = 1/2 + d/(2*l2); 
                end
                if a1<=a2 && a2<b2 && b2<=b1
                    M = 1/2 + d/(2*l1); 
                end
            else
                M = 0.5;
            end
        end
        
        % This computes the preference index for two
        % trapezoidal_fuzzy_numbers using numerical integration
        function M_APPROX =  M_approx_tfn(alpha_start,alpha_end,n,fA,fB)
            xAL = fA(1);
            nxAL =fA(2);
            nxAR = fA(3);
            xAR = fA(4);
            xBL = fB(1);
            nxBL = fB(2);
            nxBR = fB(3);
            xBR = fB(4);
            delta = (alpha_end-alpha_start)/n;
            S = 0;
            for i=1:n 
                alpha = alpha_start+i*delta; 
                a1 = alpha*nxAL+(1-alpha)*xAL;
                b1 = alpha*nxAR+(1-alpha)*xAR;
                a2 = alpha*nxBL+(1-alpha)*xBL;
                b2 = alpha*nxBR+(1-alpha)*xBR;
                Mi = pbp_class.M_interval([a1,b1],[a2,b2]);
                S = S+Mi;
            end 
            M_APPROX = S/n*(alpha_end-alpha_start);
        end
        
        function M = overlapping_1(alpha1,alpha2,fA,fB)
            % overlapping: aA < aB
            xAL = fA(1);
            nxAL =fA(2);
            nxAR = fA(3);
            xAR = fA(4);
            xBL = fB(1);
            nxBL = fB(2);
            nxBR = fB(3);
            xBR = fB(4);
            
            A0 = nxAR - xAR - (nxBL-xBL);
            B0 = xAR-xBL;
            A1 = 0.5*(nxAR-xAR-(nxAL-xAL));
            B1 = 0.5*(xAR-xAL);
            A2 = 0.5*(nxBR - xBR - (nxBL - xBL));
            B2 = 0.5*(xBR - xBL);
            if A1*B2-A2*B1 == 0
                A1 = A1+0.000001;
            end
            if A1*alpha2+B1 == 0
                A1 = A1 + 0.000001;
            end
            if A2*alpha2+B2 == 0
                A2 = A2 + 0.000001;
            end

            t1 = (1-A0^2/(8*A1*A2))*(alpha2-alpha1);
            t2 = (A0*B1 - A1*B0)^2/(8*A1^2*(A1*B2-A2*B1))*log(abs((A1*alpha2+B1)/(A1*alpha1+B1)));
            t3 = (A0*B2 - A2*B0)^2/(8*A2^2*(A1*B2-A2*B1))*log(abs((A2*alpha2+B2)/(A2*alpha1+B2)));
            M = t1-t2+t3;
        end
        
        function M = overlapping_2(alpha1,alpha2,fA,fB)
            % overlapping: aA > aB
            xAL = fA(1);
            nxAL =fA(2);
            nxAR = fA(3);
            xAR = fA(4);
            xBL = fB(1);
            nxBL = fB(2);
            nxBR = fB(3);
            xBR = fB(4);
            
            A0 = nxBR - xBR - (nxAL-xAL);
            B0 = xBR-xAL;
            A1 = 0.5*(nxAR-xAR-(nxAL-xAL));
            B1 = 0.5*(xAR-xAL);
            A2 = 0.5*(nxBR - xBR - (nxBL - xBL));
            B2 = 0.5*(xBR - xBL);
            if A1*B2-A2*B1 == 0
                A1 = A1 + 0.000001;
            end
            if A1*alpha2+B1 == 0
                A1 = A1 + 0.000001;
            end
            if A2*alpha2+B2 == 0
                A2 = A2 + 0.000001;
            end

            t1 = A0^2/(8*A1*A2)*(alpha2-alpha1);
            t2 = (A0*B1 - A1*B0)^2/(8*A1^2*(A1*B2-A2*B1))*log(abs((A1*alpha2+B1)/(A1*alpha1+B1)));
            t3 = (A0*B2 - A2*B0)^2/(8*A2^2*(A1*B2-A2*B1))*log(abs((A2*alpha2+B2)/(A2*alpha1+B2)));
            M = t1+t2-t3;
        end
        
        function M = inclusion_1(alpha1,alpha2,fA,fB)
            % A is subset of B
            xAL = fA(1);
            nxAL =fA(2);
            nxAR = fA(3);
            xAR = fA(4);
            xBL = fB(1);
            nxBL = fB(2);
            nxBR = fB(3);
            xBR = fB(4);
            
            A0 = 0.5*((nxBR - xBR) + (nxBL - xBL) - (nxAL-xAL) - (nxAR - xAR));
            B0 = 0.5*(xBL+xBR-xAL-xAR);
            A1 = 0.5*(nxBR-xBR-(nxBL-xBL));
            B1 = 0.5*(xBR-xBL);
            if A1*alpha2+B1 == 0
                A1 = A1 + 0.000001;
            end
                              
            t1 = 0.5*(alpha2-alpha1)*(1+A0/A1);
            t2 = (A0*B1 - A1*B0)/(2*A1^2)*log(abs((A1*alpha2+B1)/(A1*alpha1+B1)));
            M = t1-t2;
        end
        
        function M = inclusion_2(alpha1,alpha2,fA,fB)
             % B is subset of A
            xAL = fA(1);
            nxAL =fA(2);
            nxAR = fA(3);
            xAR = fA(4);
            xBL = fB(1);
            nxBL = fB(2);
            nxBR = fB(3);
            xBR = fB(4);
            
            A0 = 0.5*((nxBR - xBR) + (nxBL - xBL) - (nxAL-xAL) - (nxAR - xAR));
            B0 = 0.5*(xBL+xBR-xAL-xAR);
            A1 = 0.5*(nxAR-xAR-(nxAL-xAL));
            B1 = 0.5*(xAR-xAL);
            if A1*alpha2+B1 == 0
                A1 = A1 + 0.000001;
            end
            
            t1 = 0.5*(alpha2-alpha1)*(1+A0/A1);
            t2 = (A0*B1 - A1*B0)/(2*A1^2)*log(abs((A1*alpha2+B1)/(A1*alpha1+B1)));
            M = t1-t2;
        end
        
        function M = preceding_1(alpha1,alpha2,fA,fB)
            % B precedes A
            M = 0;
        end
        
         function M = preceding_2(alpha1,alpha2,fA,fB)
            % A precedes B
             M = alpha2-alpha1;
         end
        
         % This is the analytic computation of M
         function M = compute_M(fA,fB)
            xAL = fA(1);
            nxAL =fA(2);
            nxAR = fA(3);
            xAR = fA(4);
            xBL = fB(1);
            nxBL = fB(2);
            nxBR = fB(3);
            xBR = fB(4);
            
            alpha_vector = [0];
            a = xAL;
            b = nxAL;
            c = xBL;
            d = nxBL;
            alpha = (a-c)/(a-b+d-c);
            if  alpha <1 && alpha>0
                alpha_vector = [alpha_vector,alpha];
            end
            a = xAR;
            b = nxAR;
            c = xBR;
            d = nxBR;
            alpha = (a-c)/(a-b+d-c);
            if  alpha <1 && alpha>0
                alpha_vector = [alpha_vector,alpha];
            end
            a = xAR;
            b = nxAR;
            c = xBL;
            d = nxBL;
            alpha = (a-c)/(a-b+d-c);
            if  alpha <1 && alpha>0
                alpha_vector = [alpha_vector,alpha];
            end
            a = xAL;
            b = nxAL;
            c = xBR;
            d = nxBR;
            alpha = (a-c)/(a-b+d-c);
            if  alpha <1 && alpha>0
                alpha_vector = [alpha_vector,alpha];
            end
            alpha_vector = [alpha_vector,1];
            n = size(alpha_vector,2);
            alpha_star = zeros(1,n-1);
            M=0;
            for i=1:n-1 
                alpha_star(i) = (alpha_vector(i) + alpha_vector(i+1))/2;
                aA = alpha_star(i)*nxAL+(1-alpha_star(i))*xAL;
                bA = alpha_star(i)*nxAR+(1-alpha_star(i))*xAR;
                aB = alpha_star(i)*nxBL+(1-alpha_star(i))*xBL;
                bB = alpha_star(i)*nxBR+(1-alpha_star(i))*xBR;
                if aB<bB && bB<aA && aA<bA
                    M = M + pbp_class.preceding_1(alpha_vector(i),alpha_vector(i+1), fA,fB);
                end
                if aA<bA && bA<aB && aB<bB
                    M = M + pbp_class.preceding_2(alpha_vector(i),alpha_vector(i+1), fA,fB);
                end
                if aA<aB && aB<bA && bA<bB
                    M = M + pbp_class.overlapping_1(alpha_vector(i),alpha_vector(i+1), fA,fB);
                end
                if aB<aA && aA<bB && bB<bA
                    M = M + pbp_class.overlapping_2(alpha_vector(i),alpha_vector(i+1), fA,fB);
                end
                if aB<=aA && aA<bA && bA<=bB
                    M = M + pbp_class.inclusion_1(alpha_vector(i),alpha_vector(i+1), fA,fB);
                elseif aA<=aB && aB<bB && bB<=bA
                    M = M + pbp_class.inclusion_2(alpha_vector(i),alpha_vector(i+1), fA,fB);
                end
            end     
          
         end
   
    end
end

