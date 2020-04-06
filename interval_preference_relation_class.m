classdef interval_preference_relation_class < handle
    
       
    methods(Static)    
        
        % This computes the preference index for two intervals 
        function M = compute_preference_index(I_1,I_2)
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
        function M_APPROX =  compute_preference_index_for_tfns(alpha_start,alpha_end,n,fA,fB)
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
            
            for alpha=alpha_start:delta:alpha_end
                a1 = alpha*nxAL+(1-alpha)*xAL;
                b1 = alpha*nxAR+(1-alpha)*xAR;
                a2 = alpha*nxBL+(1-alpha)*xBL;
                b2 = alpha*nxBR+(1-alpha)*xBR;
                Mi = interval_preference_relation_class.compute_preference_index([a1,b1],[a2,b2]);
                S = S+Mi;
            end
            
            M_APPROX = S/n*(alpha_end-alpha_start);
        end  
        
        % This computes the preference intensity index matrix
         function M_matrix = get_preference_intensity_index_matrix_fuzzy_numbers(fuzzy_numbers)
            m = size(fuzzy_numbers,1);
            M_matrix = zeros(m,m);
            for i=1:m
                for j=i:m
                    M_matrix(i,j) = min([interval_preference_relation_class.compute_preference_index_for_tfns(0,1,10000,fuzzy_numbers(i,:),fuzzy_numbers(j,:)),1]);
                    M_matrix(j,i) = 1-M_matrix(i,j);
                end
            end
            for i=1:m
                M_matrix(i,i) = 0.5;
            end
         end
        
        % This orders the input fuzzy numbers for a given delta
        function [ordered_fuzzy_numbers,ordered_fuzzy_number_IDs] = order_fuzzy_numbers(delta,fuzzy_numbers)
            fuzzy_number_IDs = [1:size(fuzzy_numbers),1]';
            n = size(fuzzy_numbers,1);
            PM = interval_preference_relation_class.get_preference_intensity_index_matrix_fuzzy_numbers(fuzzy_numbers)
            n_i  = sum(PM>0.5+delta,2);
            p_i = n-n_i;
            n_i = n_i'
            p_i = p_i'
            k = 1:size(unique(p_i),2)
            u_k = unique(p_i)
            o_i = zeros(1,size(p_i,2));
            for i=1:size(p_i,2)
                o_i(1,i) = find(u_k==p_i(1,i));
            end
            o_i
            [~,order_index] = sort(p_i);
            %order_index
            ordered_fuzzy_numbers = fuzzy_numbers(order_index,:);
            ordered_fuzzy_number_IDs = fuzzy_number_IDs(order_index);
        end        
        
        % This plots the fuzzy numbers
        function plot_fuzzy_numbers(fuzzy_numbers,fuzzy_IDS,file_name,type)
            fig = figure();
            k = 50;
            l = 0;
            if strcmp(type,'ordered')
                hold on;
                rectangle('Position',[fuzzy_numbers(5,1) 2*(k+10) fuzzy_numbers(5,4) - fuzzy_numbers(5,1) 3*(k+10)],'linewidth',1,'edgecolor',[0 0 0],...
                    'facecolor',[0.9 0.9 0.9],'curvature',0.2,'linestyle','--');
                rectangle('Position',[fuzzy_numbers(7,1) 5*(k+10) fuzzy_numbers(6,4) - fuzzy_numbers(7,1) 2*(k+10)],'linewidth',1,'edgecolor',[0 0 0],...
                    'facecolor',[0.9 0.9 0.9],'curvature',0.2,'linestyle','--');
            end
            hold on
            numOfFuzzyNumbers = size(fuzzy_numbers,1);
            drawcolor = [0.6 0.6 0.6];
            for i = 1:numOfFuzzyNumbers
                plot(fuzzy_numbers(i,:),[l,l+k,l+k,l],'color',drawcolor,'LineWidth',3.0);
                line([0,max(max(fuzzy_numbers))+10],[l,l],'color',[0 0 0 ], 'LineWidth',0.5);
                text((fuzzy_numbers(i,1)+fuzzy_numbers(i,4))/2,l+k/2,strcat('$A_{',num2str(fuzzy_IDS(i)),'}$'),'interpreter','LaTeX','fontsize',14);
                axis([0,max(max(fuzzy_numbers))+10,0,numOfFuzzyNumbers*(k+10)])
                l = l+k+10;
            end
            hold off
            grid on
            set(gca,'ytick',[]);
            xlabel('$x$','interpreter','LaTeX','fontsize',14);
            ylabel('$\mu(x)$','interpreter','LaTeX','fontsize',14);
            
            print(fig,file_name,'-depsc','-painters');
        end
         
          % This function identifies the delta for wich the relation is transitive
          function delta = get_intransitive_triples(fuzzy_sets)
             n = size(fuzzy_sets,1);
             combin = nchoosek(1:n,3);
             delta = 0.01;
             found_intrans_triples = true;
             while found_intrans_triples
                found_intrans_triples = false; 
                delta = delta + 0.01;
                intrans_triples = [];
                 for i=1:size(combin,1)
                    permut = perms(combin(i,:));
                    for j=1:size(permut,1)
                        fA = fuzzy_sets(permut(j,1),:);
                        fB = fuzzy_sets(permut(j,2),:);
                        fC = fuzzy_sets(permut(j,3),:);
                        M_AB = interval_preference_relation_class.compute_preference_index_for_tfns(0,1,100,fA,fB);
                        M_BC = interval_preference_relation_class.compute_preference_index_for_tfns(0,1,100,fB,fC);
                        M_AC = interval_preference_relation_class.compute_preference_index_for_tfns(0,1,100,fA,fC);
                        if M_AB>=0.5+delta && M_BC>=0.5+delta && M_AC<0.5+delta
                            intrans_triples = [intrans_triples; permut(j,:)]
                            found_intrans_triples = true;
                        end    
                   end
                 end
             end
          end
 
    end
end

