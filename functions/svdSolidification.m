function [R,t,T]=svdSolidification(p_src,p_tgt)
% this is the paged version of svdm
            m_src=mean(p_src,1);
            m_tgt=mean(p_tgt,1);
            p_src=p_src-m_src; % p_src centrato in zero
            p_tgt=p_tgt-m_tgt; % p_tgt centrato in zero
            C=pagemtimes(p_src,'transpose',p_tgt,'none');
            
            nf=size(C,3);
            for i=nf:-1:1
            [U,~,V]=svd(C(:,:,i));
            R(:,:,i)=(U*diag([1 1 det(U*V')])*V')';
            t(i,:)=m_tgt(:,:,i)-m_src*R(:,:,i)';
            end
            T=R;
            T(1:3,4,:)=t';
            T(4,4,:)=1;
        end
    