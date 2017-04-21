%%% params
% GLOBAL
% .scales : size filter - window
% MORPHOLOGICAL
% .strel  : form strel
% GLCM
% .averagedir    : boolean
%     if .averagedir == 1 returns average on the 4 directions,
%     else: return all the GLCM texts separate [4 dirs (0, 45, 90, 135)]
% .dist: length lag to compute GLCM: .dist < .scales !
%     default, .dist = ceil(1/4*.scales);
% .filters: name filter:
%             OC    : Openening and closing
%             TH    : Top hat opening and closing
%             OCR   : Opening and closing by reconstruction
%             OCRTH : Top hat on OCR
%             TXT   : Occurrence texture description:
%                         - mean
%                         - variance
%                         - range
%                         - entropy
%             GLCM  : Co-occurrence (GLCM):
%                         - contrast (Variance, Inertia)
%                         - correlation
%                         - Angular second moment (Energy)
%                         - Homogeneity
% % EXA:
% params.filters ='OC';
% params.scales  = [2 4 6 8 10 11 15 67];
% params.formel  = 'disk';
% params.scales = [20];
% Sz = size(t1);
% I = t1;%
function I = contextualfeatures(I,params)

Sz = size(I);
%fprintf([params.filters '\n'])

if ~isfield(params,'formel')
    params.formel = 'disk';
end

if ~isfield(params,'dist')
    params.dist = floor((1/4).*params.scales);
end

if ~isfield(params,'averagedir')
    params.averagedir = 1;
end

if params.averagedir ~= 1
    params.averagedir = 0;
end

if size(I,3) > 1
    [tt,Ipcav]=princomp(reshape(I,Sz(1)*Sz(2),size(I,3)));
    Ipca = reshape(Ipcav(:,1),Sz(1),Sz(2));
    clear tt
else
    Ipca = reshape(I,Sz(1),Sz(2));
end

% 16 bit conversion!
if ~strcmp(class(Ipca),'uint16')
    mx = max(max(Ipca));
    mi = min(min(Ipca));
    Ipca = uint16(((Ipca - mi)./ (mx - mi)).*2^16-1);
end
I = [];
SE = [];


switch params.formel
    case 'line'
        for i=1:length(params.scales)
            SE = [SE;  strel(params.formel,params.scales(i),params.angles(i))];
        end
        
    case {'disk','square','diamond'}
        
        for i=1:length(params.scales)
            
            SE = [SE;  strel(params.formel,params.scales(i))];
            
        end
        
end



    switch params.filters
        case {'band'}
            I = reshape(Ipca,Sz(1)*Sz(2),1);
            %%%%%%%%%%%%%%%%%%%%%%% MORPHOLOGICAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %-------- Opening -------------------------------------------------
        case {'op','OP'}
            I_o = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                o_temp = imopen(Ipca,SE(i));
                I_o(:,i) = reshape(o_temp,Sz(1)*Sz(2),1);
            end
            I = [I_o];
            %-------- Closing -------------------------------------------------
        case {'cl','CL'}
            I_c = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                c_temp = imclose(Ipca,SE(i));
                I_c(:,i) = reshape(c_temp,Sz(1)*Sz(2),1);
            end
            I = [I_c];
            %-------- Top hat opening -----------------------------------------
        case {'th_op','TH_OP'}
            I_o = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                o_temp = imopen(Ipca,SE(i));
                I_o(:,i) = reshape(Ipca-o_temp,Sz(1)*Sz(2),1);
            end
            I = [I_o];
            %-------- Top hat opening -----------------------------------------
        case {'th_cl','TH_CL'}
            I_c = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                c_temp = imclose(Ipca,SE(i));
                I_c(:,i) = reshape(c_temp-Ipca,Sz(1)*Sz(2),1);
            end
            I = [I_c];
            %-------- Opening by reconstruction -------------------------------
        case {'or','OR'}
            I_o = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                m_o = imerode(Ipca,SE(i));
                o_temp = imreconstruct(m_o,Ipca);
                I_o(:,i) = reshape(o_temp,Sz(1)*Sz(2),1);
            end
            I = [I_o];
            %-------- Closing by reconstruction ------------------------------
        case {'cr','CR'}
            I_c = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                m_c = imdilate(Ipca,SE(i));
                c_temp = imcomplement(imreconstruct(imcomplement(m_c), imcomplement(Ipca)));
                I_c(:,i) = reshape(c_temp,Sz(1)*Sz(2),1);
            end
            I = [I_c];
            %-------- Top hat with opening by reconstruction ------------------
        case {'th_or','TH_OR'}
            I_o = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                m_o = imerode(Ipca,SE(i));
                o_temp = imreconstruct(m_o,Ipca);
                I_o(:,i) = reshape(Ipca - o_temp,Sz(1)*Sz(2),1);
            end
            I = [I_o];
            %-------- Top hat with opening by reconstruction ------------------
        case {'th_cr','TH_CR'}
            I_c = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                m_c = imdilate(Ipca,SE(i));
                c_temp = imcomplement(imreconstruct(imcomplement(m_c), imcomplement(Ipca)));
                I_c(:,i) = reshape(c_temp - Ipca,Sz(1)*Sz(2),1);
            end
            I = [I_c];
            
            %%%%%%%%%%%%%%%%%%%%%%% ATTRIBUTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %-------- ATTRIBUTE FILTERS: area ---------------------------------% case('ATT','att')
        case {'ATT-a','att-a'} % attribute profile, area
            % ATTRIBUTE_PROFILE returns, for 3 scales, [C3 C2 C1 IMG O1 O2 O3]
            % 'a' attribute proportional of area of objects (in pixels)
            res = attribute_profile(Ipca,'a', params.scales);
            res = reshape(double(res),Sz(1)*Sz(2),size(res,3));
            I = res(:,[1:length(params.scales) length(params.scales)+2:end]);
            %-------- ATTRIBUTE FILTERS: diagonal of BB------------------------% case('ATT','att')
        case {'ATT-d','att-d'} %attribute profile, diagonal of bounding box
            % ATTRIBUTE_PROFILE returns, for 3 scales, [C3 C2 C1 IMG O1 O2 O3]
            % 'd' attribute proportional of size of objects (diagonal in pixels)
            res = attribute_profile(Ipca,'d', params.scales);
            res = reshape(double(res),Sz(1)*Sz(2),size(res,3));
            I = res(:,[1:length(params.scales) length(params.scales)+2:end]);
            %-------- ATTRIBUTE FILTERS: inertia ------------------------------% case('ATT','att')
        case {'ATT-i','att-i'} %attribute profile, moment of intertia
            % ATTRIBUTE_PROFILE returns, for 3 scales, [C3 C2 C1 IMG O1 O2 O3]
            % 'i' attribute \in [0.1, ..., 1]
            res = attribute_profile(Ipca,'i', params.scales);
            res = reshape(double(res),Sz(1)*Sz(2),size(res,3));
            I = res(:,[1:length(params.scales) length(params.scales)+2:end]);
            %-------- ATTRIBUTE FILTERS: standard dev -------------------------% case('ATT','att')
        case {'ATT-s','att-s'} %attribute profile, standard deviation
            % ATTRIBUTE_PROFILE returns, for 3 scales, [C3 C2 C1 IMG O1 O2 O3]
            % 's' attribute proportional of std of image ?
            res = attribute_profile(Ipca,'s', params.scales);
            res = reshape(double(res),Sz(1)*Sz(2),size(res,3));
            I = res(:,[1:length(params.scales) length(params.scales)+2:end]);
            
            %%%%%%%%%%%%%%%%%%%%%%% TEXTURE OCCURRENCE %%%%%%%%%%%%%%%%%%%%%%%%
            %-------- TXT: Local average --------------------------------------%    case {'TXT','txt'}
        case {'avg','AVG'}
            I_avg = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                res = imfilter(Ipca, fspecial('average',params.scales(i)), 'replicate');
                I_avg(:,i) = reshape(res,Sz(1)*Sz(2),1);
                clear res
            end
            I = I_avg;
            %-------- TXT: Local entropy --------------------------------------%    case {'TXT','txt'}
        case {'ent','ENT'}
            I_ent = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                res = entropyfilt(Ipca, ones(params.scales(i)));%!!! occhio che lo vuole in uint16!
                I_ent(:,i) = reshape(res,Sz(1)*Sz(2),1);
                clear res
            end
            I = I_ent;
            %-------- TXT: Local standard deviation ---------------------------%    case {'TXT','txt'}
        case {'std','STD'}
            I_std = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                res = stdfilt(Ipca, true(params.scales(i)));
                I_std(:,i) = reshape(res,Sz(1)*Sz(2),1);
                clear res
            end
            I = I_std;
            %-------- TXT: Local range ----------------------------------------%    case {'TXT','txt'}
        case {'ran','RAN'}
            I_ran = zeros(Sz(1)*Sz(2),length(params.scales));
            for i=1:length(params.scales)
                res = rangefilt(Ipca, true(params.scales(i)));
                I_ran(:,i) = reshape(res,Sz(1)*Sz(2),1);
                clear res
            end
            I = I_ran;
            %%%%%%%%%%%%%%%%%%%%%%% COMMENTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         % possibile implementarne altre con handle di funzioni dichiarate
            %         % 'online' :
            %         FUN = @(x) mean(x(:));
            %         I_filt = nlfilter(Ipca, [params.scales(i) params.scales(i)], FUN);
            %         figure; imagesc(I_filt); axis image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLCM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (NOT DEBUGGED PROPERLY)
            %-------- GLCM: Contrast ------------------------------------------% case('GLCM','glcm')
            %         case {'glcm_cont','GLCM_CONT'}
            %             params.averagedir = 1;
            %             warning('only average among directions supported. indexing in IKL!') %#ok<WNTAG>
            %             matlabpool open local
            %             for i = 1:length(params.scales)
            %                 I_cont_T = zeros(Sz(1)*Sz(2),4);
            %
            %                 DISTS = [0  params.dist(i);          % 0???
            %                     -params.dist(i)  params.dist(i); % 45???
            %                     -params.dist(i) 0;               % 90???
            %                     -params.dist(i) -params.dist(i)];% 135???
            %
            %                 Ipca_exp = [Ipca(:,floor(params.scales(i)/2):-1:1) Ipca ...
            %                     Ipca(:,end-1:-1:end-floor(params.scales(i)/2))];
            %                 Ipca_exp = [Ipca_exp(floor(params.scales(i)/2):-1:1,:); ...
            %                     Ipca_exp; Ipca_exp(end-1:-1:end-floor(params.scales(i)/2),:)];
            %
            %                 Impa = im2patch(Ipca_exp,params.scales(i),1);
            %                 %             tic
            %                 parfor jj = 1:(Sz(1)*Sz(2))
            %                     temp  = reshape(Impa(jj,:),params.scales(i),params.scales(i));
            %                     [GLCM,SI]  = graycomatrix(temp,'GrayLimits',[],...
            %                         'NumLevels',8,'Offset',DISTS,'Symmetric',true);
            %                     stats = graycoprops(GLCM,'Contrast');
            
        otherwise
            
            fprintf(['Unkonwn filter. Quitting! :s \n valid filters'])
            
            return
            
            
    end
    

I = double(I);
