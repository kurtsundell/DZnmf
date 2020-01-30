%% DZnmf %%

function varargout = DZnmf_v05(varargin)
% DZNMF_V05 MATLAB code for DZnmf_v05.fig
%      DZNMF_V05, by itself, creates a new DZNMF_V05 or raises the existing
%      singleton*.
%
%      H = DZNMF_V05 returns the handle to a new DZNMF_V05 or the handle to
%      the existing singleton*.
%
%      DZNMF_V05('CALLBACK',hObject,eventData,H,...) calls the local
%      function named CALLBACK in DZNMF_V05.M with the given input arguments.
%
%      DZNMF_V05('Property','Value',...) creates a new DZNMF_V05 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DZnmf_v05_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DZnmf_v05_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIH

% Edit the above text to modify the response to help DZnmf_v05

% Last Modified by GUIDE v2.5 03-Jun-2018 18:30:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DZnmf_v05_OpeningFcn, ...
                   'gui_OutputFcn',  @DZnmf_v05_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% set paths for writing spreadsheets
javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');

% End initialization code - DO NOT EDIT

function DZnmf_v05_OpeningFcn(hObject, eventdata, H, varargin)
H.output = hObject;
guidata(hObject, H);

function varargout = DZnmf_v05_OutputFcn(hObject, eventdata, H) 
varargout{1} = H.output;

set(H.sigma1,'Enable','on')  %disable checkbox
set(H.sigma2,'Enable','on')  %disable checkbox
set(H.sigma1,'Value',1);

set(H.check_nmf,'Value',1)
set(H.check_nmf_opt,'Value',0);
set(H.num_sources_opt,'Enable','off');
set(H.max_num,'Enable','off');
set(H.opt_num,'Enable','off');
set(H.opt_num_result,'Enable','off');
set(H.num_sources,'Enable','on');
set(H.srcs,'Enable','on');
set(H.cancel,'Enable','off');
set(H.replot_res,'Enable','off');

set(H.pdps,'Value',1)
set(H.kern_opt,'Enable','off');
set(H.kern_set,'Enable','off');
set(H.kern,'Enable','off');
set(H.Myr,'Enable','off');

axes(H.stacked_sink_samples);
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

axes(H.stacked_reconstructed_sources);
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])

optimized = 0;
data = [];
N = [];
name = [];
text = [];
pdps_active = [];

global cancel
cancel = 0;

H.optimized = optimized;
H.data = data;
H.N = N;
H.name = name;
H.text = text;
H.pdps_active = pdps_active;

guidata(hObject, H);

%% --- Executes during object creation, after setting all properties.
function Loaded_samples_CreateFcn(hObject, eventdata, H)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% CHECKBOX Load Ages %%
function load_ages_Callback(hObject, eventdata, H)

if get(H.load_ages,'Value') == 1
	set(H.load_dists,'Value',0);
	set(H.sigma1,'Enable','on')  %disable checkbox
	set(H.sigma2,'Enable','on')  %disable checkbox
	set(H.sigma1,'Value',1);
	set(H.xaxis_min,'Enable','on');
	set(H.xaxis_max,'Enable','on');
	set(H.xaxis_int,'Enable','on');
	set(H.minx,'Enable','on');
	set(H.maxx,'Enable','on');
	set(H.intx,'Enable','on');
	set(H.pdps,'Enable','on');
	set(H.kdes,'Enable','on');
end

%% CHECKBOX Load Distributions %%
function load_dists_Callback(hObject, eventdata, H)

if get(H.load_dists,'Value') == 1
	set(H.load_ages,'Value',0);
	set(H.sigma1,'Enable','off');
	set(H.sigma2,'Enable','off');
	set(H.sigma1,'Value',0);
	set(H.sigma2,'Value',0);
	set(H.xaxis_min,'Enable','off');
	set(H.xaxis_max,'Enable','off');
	set(H.xaxis_int,'Enable','off');
	set(H.minx,'Enable','off');
	set(H.maxx,'Enable','off');
	set(H.intx,'Enable','off');
	set(H.pdps,'Enable','off');
	set(H.kdes,'Enable','off');
end

%% CHECKBOX Run NMF %%
function check_nmf_Callback(hObject, eventdata, H)

if get(H.check_nmf,'Value') == 1
	set(H.check_nmf_opt,'Value',0);
	set(H.num_sources_opt,'Enable','off');
	set(H.max_num,'Enable','off');
	set(H.opt_num,'Enable','off');
	set(H.opt_num_result,'Enable','off');
	set(H.num_sources,'Enable','on');
	set(H.srcs,'Enable','on');
	set(H.cancel,'Enable','off');
	set(H.replot_res,'Enable','off');
end


%% CHECKBOX Run NMF and find optimal number of sources %%
function check_nmf_opt_Callback(hObject, eventdata, H)

if get(H.check_nmf_opt,'Value') == 1
	set(H.check_nmf,'Value',0);
	set(H.num_sources_opt,'Enable','on');
	set(H.max_num,'Enable','on');
	set(H.opt_num,'Enable','on');
	set(H.opt_num_result,'Enable','on');
	set(H.num_sources,'Enable','off');
	set(H.srcs,'Enable','off');
	set(H.cancel,'Enable','on');
	set(H.replot_res,'Enable','on');
end

%% CHECKBOX PDPs %%
function pdps_Callback(hObject, eventdata, H)

if get(H.pdps,'Value') == 1
	set(H.kdes,'Value',0)
	set(H.kern_opt,'Enable','off');
	set(H.kern_opt,'Value',0);
	set(H.kern_set,'Value',0);
	set(H.kern_set,'Enable','off');
	set(H.kern,'Enable','off');
	set(H.Myr,'Enable','off');
end
%{
cla(H.stacked_reconstructed_sources,'reset');
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

pdps_active = [];
name = H.name;

set(H.Active_samples, 'String', []);
set(H.Loaded_samples, 'String', name);
set(H.Active_samples, 'Value',1);

ylabel([])
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])
%}
%active_selected = {1};
%loaded_selected = {1};
%H.active_selected = active_selected;
%H.loaded_selected = loaded_selected;
%H.pdps_active = pdps_active;
guidata(hObject, H);

%% CHECKBOX KDEs %%
function kdes_Callback(hObject, eventdata, H)

if get(H.kdes,'Value') == 1
	set(H.pdps,'Value',0)
	set(H.kern_opt,'Enable','on');
	set(H.kern_opt,'Value',1);
	set(H.kern_set,'Value',0);
	set(H.kern_set,'Enable','on');
	set(H.kern,'Enable','on');
end
%{
cla(H.stacked_reconstructed_sources,'reset');
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

pdps_active = [];
name = H.name;

set(H.Active_samples, 'String', []);
set(H.Loaded_samples, 'String', name);
set(H.Active_samples, 'Value',1);

ylabel([])
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])
%}
%active_selected = {1};
%loaded_selected = {1};
%H.active_selected = active_selected;
%H.loaded_selected = loaded_selected;
%H.pdps_active = pdps_active;
guidata(hObject, H);

%% PUSHBUTTON Reset to Default Parameters
function reset_params_Callback(hObject, eventdata, H)

set(H.iterations,'String','10000')
set(H.Final_Residual,'String','1E-8');
set(H.tof,'String','2E-16');

%% PUSHBUTTON Load Data %%
function Load_Callback(hObject, eventdata, H)

pdps_active = H.pdps_active;

if isempty(pdps_active) == 0
	
cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

cla(H.stacked_reconstructed_sources,'reset');
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])

pdps_active = [];
name = H.name;

set(H.Active_samples, 'String', []);
set(H.Loaded_samples, 'String', name);
set(H.Active_samples, 'Value',1);

ylabel([])
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])

end




data = H.data;
N = H.N;
name = H.name;
text = H.text;

if isempty(data) == 0
	[filename, pathname] = uigetfile({'*'},'File Selector');
	fullpathname = strcat(pathname, filename);
	[numbers, text_new, data_new_tmp] = xlsread(fullpathname);
	if isempty(text_new) == 1
		err_dlg=errordlg('Please add headers (names) to your input data for both ages and uncertainties. Press Example Input button to see an example.','Error!');
		waitfor(err_dlg);
	else
		new_path = [get(H.filepath,'String'), '   &   ', fullpathname];
		set(H.filepath,'String',new_path);
		data_new = cell2mat(data_new_tmp(2:end,:));
		data_new_length = max(length(data_new(:,1)),length(data(:,1)));
		data_new_N = length(data_new(1,:))/2;
		data_tmp = data;
		data = zeros(data_new_length,N*2+data_new_N*2);
		data(1:length(data_tmp(:,1)),1:N*2) = data_tmp;
		data(1:length(data_new(:,1)),N*2+1:N*2+data_new_N*2) = data_new;
		N = N + data_new_N;
		for i = 1:data_new_N
			name_new(i,1) = text_new(1,i*2-1);
		end
		name = [name;name_new];
		text = [text,text_new];
	end
end

if isempty(data) == 1
	[filename, pathname] = uigetfile({'*'},'File Selector');
	fullpathname = strcat(pathname, filename);
	[numbers, text, data_tmp] = xlsread(fullpathname);
	if isempty(text) == 1
		err_dlg=errordlg('Please add headers (names) to your input data for both ages and uncertainties. Press Example Input button to see an example.','Error!');
		waitfor(err_dlg);
	else
		set(H.filepath, 'String', fullpathname);
		data = cell2mat(data_tmp(2:end,:));
		[dataR,dataC]=size(data);
		N = (dataC/2); % Number of source samples
		clear name
	end
end
	

	
	for i = 1:N
		name(i,1) = text(1,i*2-1);
	end




rad_on=get(H.input_data,'selectedobject');
switch rad_on
	case H.sigma2
	for i = 1:N
		data(:,i*2) = data(:,i*2)./2;
	end
end

set(H.Loaded_samples, 'String', name);

pdps_active = [];
set(H.Loaded_samples,'Value', 1);

active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
H.name = name;
H.text = text;
H.N = N;
H.data = data;

guidata(hObject, H);

%% PUSHBUTTON Run NMF! %%
function run_nmf_Callback(hObject, eventdata, H)
opengl hardware
cla(H.stacked_reconstructed_sources,'reset');
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])
%cla(H.Final_Residual_v_Rank,'reset');
%cla(H.SSR_v_Rank,'reset');

c = cellstr(get(H.Active_samples,'String'));
check = c{1,1};

if isempty(check) == 1
	err_dlg=errordlg('There are no active samples!','But wait...');
	waitfor(err_dlg);
else

	set(H.xaxis_min,'Enable','off');
	set(H.xaxis_max,'Enable','off');
	set(H.xaxis_int,'Enable','off');
	set(H.pdps,'Enable','off');
	set(H.kdes,'Enable','off');
	
pdps_active = H.pdps_active;

name = H.name;

option.iterations = str2num(get(H.iterations,'String'));
option.residual = str2num(get(H.Final_Residual,'String'));
option.tof = str2num(get(H.tof,'String'));

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));

x = x_min:x_int:x_max;

global iter_num
global cancel

cancel = 0;

if get(H.check_nmf, 'Value') == 1
	nsources = str2num(get(H.num_sources,'String'));
	iter_num = nsources;
	[source_PDP,Weightings,numIter,tElapsed,finalResidual]=DZnmf(pdps_active,nsources,option);
    
    %%Saylor edits 30 May
    while numIter==str2num(get(H.iterations,'String'))
        answer=questdlg('Factorization did not converge before maximum number of iterations was reached.  Run multiple iterations?', ...
            'Maximum number of iterations reached','Yes','No','Yes');
        switch answer
            case 'Yes'
               prompt='Input number of iterations';
               title='Input';
               dims=[1 50];
               default_input={'50'};
               repeat_iterations=str2num(cell2mat(inputdlg(prompt,title,dims,default_input)));
               for i=1:repeat_iterations
                   [source_PDP(:,:,i),Weightings(:,:,i),numIter(i),tElapsed(i),finalResidual(i)]=DZnmf(pdps_active,nsources,option);
               end
               [minimum_finalResidual, min_index]=min(finalResidual);
               source_PDP=source_PDP(:,:,min_index);
               Weightings=Weightings(:,:,min_index);
               numIter=numIter(min_index);
               finalResidual=finalResidual(min_index);
            case 'No'
                opts.WindowStyle='modal';
                opts.Interpreter = 'none';
                f=errordlg('Treat results with caution', 'Error',opts);
                break
        end
    end
    %%END Saylor edits 30 May
    
    %%Saylor edits 29 May
    %reconstruct sink age distributions and CDFs
    reconstructed_sinks_PDP=source_PDP*Weightings;
    reconstructed_sinks_CDF=cumsum(reconstructed_sinks_PDP);
    cdfs_active=cumsum(pdps_active);
    %run metrics comparing input to reconstructed sinks
    %Cross-correlation
    [trash,number_of_sinks]=size(reconstructed_sinks_PDP);
    for i=1:number_of_sinks
        R2(1,i)=r2(reconstructed_sinks_PDP(:,i), pdps_active(:,i));
    end
    mean_R2=mean(R2);
    stdev_R2=std(R2);
    set(H.mean_R2, 'String', mean_R2);
    set(H.stdev_R2, 'String', stdev_R2);
    
    %Likeness
    for i=1:number_of_sinks
        Likeness(1,i)=like(reconstructed_sinks_PDP(:,i), pdps_active(:,i));
    end
    mean_Likeness=mean(Likeness);
    stdev_Likeness=std(Likeness);
    set(H.mean_Likeness, 'String', mean_Likeness);
    set(H.stdev_Likeness, 'String', stdev_Likeness);
    
    %Kuiper
    for i=1:number_of_sinks
        Kuiper(1,i)=max(reconstructed_sinks_CDF(:,i)-cdfs_active(:,i))+max(cdfs_active(:,i)-reconstructed_sinks_CDF(:,i));
    end
    mean_Kuiper=mean(Kuiper);
    stdev_Kuiper=std(Kuiper);
    set(H.mean_Kuiper, 'String', mean_Kuiper);
    set(H.stdev_Kuiper, 'String', stdev_Kuiper);
    
    %KS
    for i=1:number_of_sinks
        K_S(1,i)=max(abs(reconstructed_sinks_CDF(:,i)-cdfs_active(:,i)));
    end
    mean_K_S=mean(K_S);
    stdev_K_S=std(K_S);
    set(H.mean_K_S, 'String', mean_K_S);
    set(H.stdev_K_S, 'String', stdev_K_S);
    
    %Final Residual
    set(H.mean_finalResidual, 'String', finalResidual);
    %%END Saylor edits 29 May
    
end

if get(H.check_nmf_opt, 'Value') == 1		
	nsources = str2num(get(H.num_sources_opt,'String'));
	[source_PDP_loop,Weightings_loop,finalResidual_loop, coefficient_count_loop,R2_loop, numIter_loop]=DZnmf_loop(pdps_active,nsources);
	xdata=transpose(2:1:nsources);
	ydata=transpose(finalResidual_loop);
        
	%test integer breakpoint
	i=1;
	clear breakpoint brkpnt_coef SSR
	for xb=2:max(xdata)
		P0=[2 2 2];
		plusfun = @(x,xb)min(x,xb);
		model = @(P,x) P(1) + P(2)*x + P(3)*(plusfun(x,xb)-xb);
		Pfit = lsqcurvefit(model,P0,xdata,ydata);
		modelpred = model(Pfit,sort(xdata));
		squared_residuals=(modelpred-ydata).^2;
		SSR(i)=max(cumsum(squared_residuals));
		i=i+1;
	end
	[trash,brkpnt_coef]=min(SSR);
	breakpoint=xdata(brkpnt_coef);

    
    %%Saylor edits 29 May
    %reconstruct sink age distributions and CDFs
    for i=1:nsources-1
        reconstructed_sinks_PDP(:,:,i)=source_PDP_loop(:,:,i)*Weightings_loop(:,:,i);
        reconstructed_sinks_CDF(:,:,i)=cumsum(reconstructed_sinks_PDP(:,:,i));
    end
    cdfs_active=cumsum(pdps_active);
    
    %run metrics comparing input to reconstructed sinks
    [trash,number_of_sinks, number_of_loops]=size(reconstructed_sinks_PDP);
    %Cross-correlation
    for j=1:number_of_loops
        for i=1:number_of_sinks
            R2(j,i)=r2(reconstructed_sinks_PDP(:,i,j), pdps_active(:,i));
        end
        mean_R2(j)=mean(R2(j,:));
        stdev_R2(j)=std(R2(j,:));
    end
    
    %Likeness
    for j=1:number_of_loops
        for i=1:number_of_sinks
            Likeness(j,i)=like(reconstructed_sinks_PDP(:,i,j), pdps_active(:,i));
        end
        mean_Likeness(j)=mean(Likeness(j,:));
        stdev_Likeness(j)=std(Likeness(j,:));
    end
    %Kuiper
    for j=1:number_of_loops
        for i=1:number_of_sinks
            Kuiper(j,i)=max(reconstructed_sinks_CDF(:,i,j)-cdfs_active(:,i))+max(cdfs_active(:,i)-reconstructed_sinks_CDF(:,i,j));
        end
        mean_Kuiper(j)=mean(Kuiper(j,:));
        stdev_Kuiper(j)=std(Kuiper(j, :));
    end
    %KS
    for j=1:number_of_loops
        for i=1:number_of_sinks
        K_S(j,i)=max(abs(reconstructed_sinks_CDF(:,i,j)-cdfs_active(:,i)));
        end
        mean_K_S(j)=mean(K_S(j,:));
        stdev_K_S(j)=std(K_S(j,:));
    end
    %%END Saylor edits 29 May
        
    source_PDP = source_PDP_loop(:,:,breakpoint);
	nsources = breakpoint;
    
    %%Saylor edits 29 May
    set(H.mean_R2, 'String', mean_R2(breakpoint));
    set(H.stdev_R2, 'String', stdev_R2(breakpoint));
    set(H.mean_Likeness, 'String', mean_Likeness(breakpoint));
    set(H.stdev_Likeness, 'String', stdev_Likeness(breakpoint));
    set(H.mean_Kuiper, 'String', mean_Kuiper(breakpoint));
    set(H.stdev_Kuiper, 'String', stdev_Kuiper(breakpoint));
    set(H.mean_K_S, 'String', mean_K_S(breakpoint));
    set(H.stdev_K_S, 'String', stdev_K_S(breakpoint));
    set(H.mean_finalResidual, 'String', finalResidual_loop(breakpoint));
    %%END Saylor edits 29 May
    
end



axes(H.stacked_reconstructed_sources);

if nsources == 1
	plot(x,source_PDP, 'k','linewidth',1.5);
else
	for i = 1:nsources
		dist_max(i,1) = max(source_PDP(:,i));
	end
	for i = 1:nsources
		dist_max2(i,1) = sum(dist_max(1:i,1));
	end
	pdp_sum = sum(dist_max2);
	for i = 1:nsources
		pdp_adj(:,i) = source_PDP(:,i) - dist_max2(i,1) + dist_max2(end,1);
	end
	for i = 1:nsources
		plot(x,pdp_adj(:,i), 'k','linewidth',1.5);
	hold on
	end
	axis([x_min x_max 0 dist_max2(end,1)])
end

xlabel('Age (Ma)')

ylabel('Relative Probability')
ytickformat('%.0f')
set(H.stacked_reconstructed_sources,'ytick',[]);


if get(H.check_nmf_opt, 'Value') == 1	
	figure;
	xb=breakpoint;
	set(H.opt_num_result,'String',xb);
	plusfun = @(x,xb)min(x,xb);
	model = @(P,x) P(1) + P(2)*x + P(3)*(plusfun(x,xb)-xb);
	P0=[2 2 2];
	Pfit = lsqcurvefit(model,P0,xdata,ydata);
	modelpred = model(Pfit,sort(xdata));
	plot(xdata,ydata,'o',sort(xdata),modelpred,'r-');
	%title('Final Residual versus Rank')
	xlabel('Rank')
	ylabel('Final Residual')

	%plot SSR vs rank
	%axes(H.SSR_v_Rank);
	figure;
	plot(xdata,SSR,'-o');
	%concatenate data for output
	output=horzcat(xdata,ydata,modelpred,transpose(SSR));
	%title('SSR versus Rank')
	xlabel('Rank')
	ylabel('SSR1 + SSR2')
	
	optimized = 1;
	H.optimized = optimized;
	H.xdata = xdata;
	H.ydata = ydata;
	H.SSR = SSR;
	H.modelpred = modelpred;
end

if get(H.check_nmf, 'Value') == 1 
	export{length(x)+2,nsources+9+number_of_sinks} = [];
	export(1,1) = {'Reconstructed Source Distributions'};
	export(2,1) = {'Age (Ma)'};
	export(3:end,1) = num2cell(x');
	BaseName='Source ';
	for i=1:nsources
		export{2,i+1}=[BaseName,num2str(i)];
	end
	export(3:length(x)+2,2:nsources+1) = num2cell(source_PDP);
	export(1,nsources+4) = {'Compare Input & Reconstructed Sink Samples'};
	export(3,nsources+4) = {'Cross-correlation'};
	export(4,nsources+4) = {'Likeness'};
	export(5,nsources+4) = {'Kuiper V Value'};
	export(6,nsources+4) = {'KS D Value'};
	export(8,nsources+4) = {'Number of sink samples'};
    export(9,nsources+4) = {'Mean Cross-correlation of sink samples'};
    export(10,nsources+4)= {'Range Cross-correlation of sink samples'};
    export(11,nsources+4)= {'Mean Kuiper V value of sink samples'};
    export(12,nsources+4)= {'Range Kuiper V value of sink samples'};
    export(13,nsources+4)= {'Final Residual'};
    export(14,nsources+4) = {'Number of Iterations '};
	export(2,nsources+5) = {'Mean'};
	export(3,nsources+5) = num2cell(mean_R2);
	export(4,nsources+5) = num2cell(mean_Likeness);
	export(5,nsources+5) = num2cell(mean_Kuiper);
	export(6,nsources+5) = num2cell(mean_K_S);
    export(8,nsources+5) = num2cell(str2num(get(H.nsamples,'String')));
    export(9,nsources+5) = num2cell(str2num(get(H.mean_R2_sink,'String')));
    export(10,nsources+5)= num2cell(str2num(get(H.range_R2_sink,'String')));
    export(11,nsources+5)= num2cell(str2num(get(H.mean_V_sink,'String')));
    export(12,nsources+5)= num2cell(str2num(get(H.range_V_sink,'String')));
	export(13,nsources+5) = num2cell(finalResidual);
    export(14,nsources+5) = num2cell(numIter);
	export(3:6,nsources+6) = {'±'};
	export(2,nsources+7) = {'St. Dev.'};
	export(3,nsources+7) = num2cell(stdev_R2);
	export(4,nsources+7) = num2cell(stdev_Likeness);
	export(5,nsources+7) = num2cell(stdev_Kuiper);
	export(6,nsources+7) = num2cell(stdev_K_S);
	BaseName='Weights Source ';
	for i=1:nsources
		export{i+2,nsources+9}=[BaseName,num2str(i)];
	end
	export(nsources+4,nsources+9) = {'Cross-correlation'};
	export(nsources+5,nsources+9) = {'Likeness'};
	export(nsources+6,nsources+9) = {'Kuiper V Value'};
	export(nsources+7,nsources+9) = {'KS D Value'};
	export(1,nsources+10) = {'Sink Sample Weights'};
	export(3:nsources+2,nsources+10:nsources+number_of_sinks+9) = num2cell(Weightings);
	export(nsources+4,nsources+10:nsources+number_of_sinks+9) = num2cell(R2);
	export(nsources+5,nsources+10:nsources+number_of_sinks+9) = num2cell(Likeness);
	export(nsources+6,nsources+10:nsources+number_of_sinks+9) = num2cell(Kuiper);
	export(nsources+7,nsources+10:nsources+number_of_sinks+9) = num2cell(K_S);
	export(2,nsources+10:nsources+number_of_sinks+9) = cellstr(get(H.Active_samples,'String'));
end

if get(H.check_nmf_opt, 'Value') == 1
	tst{1,iter_num} = [];
	tst2{1,iter_num} = [];
	BaseName='Source ';
	BaseName2='Weights Source ';
	for i=1:iter_num
		tst{1,i}=[BaseName,num2str(i)];
		tst2{1,i}=[BaseName2,num2str(i)];
	end
	export{length(x)+2,nsources+9+number_of_sinks,iter_num-1} = [];

	for k = 1:iter_num-1
		export(1,1,k) = {'Reconstructed Source Distributions'};
		export(2,1,k) = {'Age (Ma)'};
		export(3:end,1,k) = num2cell(x');
		export(2,2:k+2,k) = tst(1,1:k+1);
		export(3:length(x)+2,2:k+2,k) = num2cell(source_PDP_loop(:,1:k+1,k));
		export(1,k+5,k) = {'Compare Input & Reconstructed Sink Samples'};
		export(3,k+5,k) = {'Cross-correlation'};
		export(4,k+5,k) = {'Likeness'};
		export(5,k+5,k) = {'Kuiper V Value'};
		export(6,k+5,k) = {'KS D Value'};
		export(8,k+5,k) = {'Number of sink samples'};
        export(9,k+5,k) = {'Mean Cross-correlation of sink samples'};
        export(10,k+5,k)= {'Range Cross-correlation of sink samples'};
        export(11,k+5,k)= {'Mean Kuiper V value of sink samples'};
        export(12,k+5,k)= {'Range Kuiper V value of sink samples'};
        export(13,k+5,k)= {'Final Residual'};
        export(14,k+5,k) = {'Number of Iterations '};
		export(2,k+6,k) = {'Mean'};
		export(3,k+6,k) = num2cell(mean_R2(1,k));
		export(4,k+6,k) = num2cell(mean_Likeness(1,k));
		export(5,k+6,k) = num2cell(mean_Kuiper(1,k));
		export(6,k+6,k) = num2cell(mean_K_S(1,k));
        export(8,k+6,k) = num2cell(str2num(get(H.nsamples,'String')));
        export(9,k+6,k) = num2cell(str2num(get(H.mean_R2_sink,'String')));
        export(10,k+6,k)= num2cell(str2num(get(H.range_R2_sink,'String')));
        export(11,k+6,k)= num2cell(str2num(get(H.mean_V_sink,'String')));
        export(12,k+6,k)= num2cell(str2num(get(H.range_V_sink,'String')));
		export(13,k+6,k) = num2cell(finalResidual_loop(1,k));
		export(14,k+6,k) = num2cell(numIter_loop(k));
		export(3:6,k+7,k) = {'±'};
		export(2,k+8,k) = {'St. Dev.'};
		export(3,k+8,k) = num2cell(stdev_R2(1,k));
		export(4,k+8,k) = num2cell(stdev_Likeness(1,k));
		export(5,k+8,k) = num2cell(stdev_Kuiper(1,k));
		export(6,k+8,k) = num2cell(stdev_K_S(1,k));
		export(3:k+3,k+10,k) = tst2(1,1:k+1);
		export(k+5,k+10,k) = {'Cross-correlation'};
		export(k+6,k+10,k) = {'Likeness'};
		export(k+7,k+10,k) = {'Kuiper V Value'};
		export(k+8,k+10,k) = {'KS D Value'};
		export(1,k+10,k) = {'Sink Sample Weights'};
		export(3:k+3,k+11:k+10+number_of_sinks,k) = num2cell(Weightings_loop(1:k+1,:,k));
		export(k+5,k+11:k+10+number_of_sinks,k) = num2cell(R2(k,:));
		export(k+6,k+11:k+10+number_of_sinks,k) = num2cell(Likeness(k,:));
		export(k+7,k+11:k+10+number_of_sinks,k) = num2cell(Kuiper(k,:));
		export(k+8,k+11:k+10+number_of_sinks,k) = num2cell(K_S(k,:));
		export(2,k+11:k+10+number_of_sinks,k) = cellstr(get(H.Active_samples,'String'));
	end
end






H.iter_num = iter_num;
H.nsources = nsources;
H.source_PDP = source_PDP;
H.export = export;

end


guidata(hObject, H);

%% PUSHBUTTON Clear All %%
function clear_all_Callback(hObject, eventdata, H)

set(H.filepath,'String',[]);
set(H.Loaded_samples,'String',[]);
set(H.Active_samples,'String',[]);
cla(H.stacked_sink_samples,'reset');
cla(H.stacked_reconstructed_sources,'reset');
%cla(H.Final_Residual_v_Rank,'reset');
%cla(H.SSR_v_Rank,'reset');
set(H.iterations,'String','10000')
set(H.Final_Residual,'String','1E-8');
set(H.tof,'String','2E-16');
set(H.opt_num_result,'String','?');

active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
pdps_active = [];
H.pdps_active = pdps_active;

axes(H.stacked_sink_samples);
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

axes(H.stacked_reconstructed_sources);
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])

set(H.xaxis_min,'Enable','on');
set(H.xaxis_max,'Enable','on');
set(H.xaxis_int,'Enable','on');

optimized = 0;

data = [];
N = [];
name = [];
text = [];

H.optimized = optimized;
H.data = data;
H.N = N;
H.name = name;
H.text = text;
%%Saylor edits 29 May
    set(H.mean_R2, 'String', 'N/A');
    set(H.stdev_R2, 'String', 'N/A');
    set(H.mean_Likeness, 'String', 'N/A');
    set(H.stdev_Likeness, 'String', 'N/A');
    set(H.mean_Kuiper, 'String', 'N/A');
    set(H.stdev_Kuiper, 'String', 'N/A');
    set(H.mean_K_S, 'String', 'N/A');
    set(H.stdev_K_S, 'String', 'N/A');
    set(H.mean_finalResidual, 'String', 'N/A');
    %%Saylor edits 29 May
guidata(hObject, H);

%% LISTBOX Loaded Samples %%
function Loaded_samples_Callback(hObject, eventdata, H)

loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_selected = {get(H.Loaded_samples,'Value')};

H.loaded_names = loaded_names;
H.loaded_selected = loaded_selected;
guidata(hObject, H);

%% LISTBOX Active Samples %%
function Active_samples_Callback(hObject, eventdata, H)

active_names = cellstr(get(H.Active_samples,'String'));
active_selected = {get(H.Active_samples,'Value')};

H.active_names = active_names;
H.active_selected = active_selected;
guidata(hObject, H);

%% PUSHBUTTON Activate Samples
function Activate_Callback(hObject, eventdata, H)
opengl hardware
pdps_active = H.pdps_active;
text = H.text;
N = H.N;
data = H.data;
loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_selected = H.loaded_selected;

if isempty(loaded_names) == 1
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
elseif length(char(loaded_names)) == 0
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
else
	
	set(H.xaxis_min,'Enable','off');
	set(H.xaxis_max,'Enable','off');
	set(H.xaxis_int,'Enable','off');
	set(H.pdps,'Enable','off');
	set(H.kdes,'Enable','off');

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));

x = x_min:x_int:x_max;
	
active_names = cellstr(get(H.Active_samples,'String'));
active_names(all(cellfun('isempty',active_names),2),:) = [];
active_names = [active_names; loaded_names(loaded_selected{1,1},1)];
set(H.Active_samples, 'String', active_names);

contents = get(H.Loaded_samples,'String');
if length(contents)<1; return; end %if already empty, do nothing
Index=get(H.Loaded_samples,'Value');
contents(Index)=[]; %remove the item
Value=Index-1;
if Value<1; Value=1;end %take care of exception
set(H.Loaded_samples,'String',contents,'Value',1);

new_names = loaded_names(loaded_selected{1,1});

for i = 1:N
	for j = 1:length(new_names)
		if strcmp(text(1,i*2-1), new_names(j,1)) == 1
			new_data_name{:,i*2-1} = new_names(j,1);
			new_data(:,i*2-1:i*2) = data(:,i*2-1:i*2);
		end
	end
end

new_data(isnan(new_data))=0;
new_data( :, all(~new_data,1) ) = [];
new_data_name = new_data_name(~cellfun('isempty',new_data_name));

if get(H.pdps,'Value') == 1
	for i = 1:length(new_data_name)
		m = nonzeros(new_data(:,i*2-1));
		m = m(isfinite(m(:,1)),:);
		s = nonzeros(new_data(:,i*2));
		s = s(isfinite(s(:,1)),:);
		f = zeros(length(m),length(x));
		for j = 1:length(m)
			f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*x_int);
		end
		pdps_new(:,i) = ((sum(f, 1))/length(m)).';
	end
end

if get(H.kdes,'Value') == 1
	rad_on_kernel=get(H.Dist_opts,'selectedobject');
	switch rad_on_kernel
		case H.kern_opt
			xA = transpose(x);
			for i = 1:length(new_data_name)
				dist_data = nonzeros(new_data(:,i*2-1));
				[bandwidth,kdeA,xmesh1,cdf]=kde(dist_data,length(x),x_min,x_max);
				pdps_new(:,i) = transpose(interp1(xmesh1, kdeA, xA));
				clear dist_data kdeA
			end
		case H.kern_set
			for i = 1:length(new_data_name)
				dist_data = nonzeros(new_data(:,i*2-1));
				kernel = str2num(get(H.kern,'String'));
				kernel_dist_data(1:length(nonzeros(new_data(:,i*2-1))),1) = kernel;
				pdps_new(:,i) = pdp5(dist_data(:,1),kernel_dist_data,x_min,x_max,x_int); 
			end
		end
end

pdps_active = [pdps_active,pdps_new];

axes(H.stacked_sink_samples);
nsamples = length(pdps_active(1,:));
if nsamples == 1
	plot(x,pdps_active, 'k','linewidth',1.5);
	axis([x_min x_max 0 max(max(pdps_active))])
else
	for i = 1:nsamples
		dist_max(i,1) = max(pdps_active(:,i));
	end
	for i = 1:nsamples
		dist_max2(i,1) = sum(dist_max(1:i,1));
	end
	pdp_sum = sum(dist_max2);
	for i = 1:nsamples
		pdp_adj(:,i) = pdps_active(:,i) - dist_max2(i,1) + dist_max2(end,1);
	end
	for i = 1:nsamples
		plot(x,pdp_adj(:,i), 'k','linewidth',1.5);
	hold on
	end
	axis([x_min x_max 0 dist_max2(end,1)])
end

xlabel('Age (Ma)')
ylabel('Relative Probability')
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])

%Mean Kuiper
if nsamples==1;
    set(H.mean_V_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_V_sink(i,j)=max(cumsum(pdps_active(:,i))-...
                cumsum(pdps_active(:,j)))+max(cumsum(pdps_active(:,j))...
                -cumsum(pdps_active(:,i)));
        end
    end
    mean_V_sink_temp=[];
    for i=1:nsamples
        mean_V_sink_temp=vertcat(mean_V_sink_temp,diag(mean_V_sink,i));
    end
    mean_V_sink=round(mean(mean_V_sink_temp),3);
    range_V_sink=round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end
    
set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;

guidata(hObject, H);

end

%% PUSHBUTTON Deactivate Samples
function Deactivate_Callback(hObject, eventdata, H)

pdps_active = H.pdps_active;
text = H.text;
N = H.N;
data = H.data;
active_names = cellstr(get(H.Active_samples,'String'));
active_selected = H.active_selected;

if isempty(active_names) == 1
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
elseif length(char(active_names)) == 0
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
else

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));

loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_names(all(cellfun('isempty',loaded_names),2),:) = [];
loaded_names = [loaded_names; active_names(active_selected{1,1},1)];
for i = 1:length(loaded_names)
	if length((strsplit(loaded_names{i,1},','))) > 1
		loaded_names(length(loaded_names)+1:length(loaded_names)+length((strsplit(loaded_names{i,1},','))),1) = strsplit(loaded_names{i,1},',');
		loaded_names{i,1} = [];
	end
end
loaded_names(all(cellfun('isempty',loaded_names),2),:) = [];

set(H.Loaded_samples, 'String', loaded_names);

contents = get(H.Active_samples,'String');
if length(contents)<1; return; end %if already empty, do nothing
Index=get(H.Active_samples,'Value');
contents(Index)=[]; %remove the item
Value=Index-1;
if Value<1; Value=1;end %take care of exception
set(H.Active_samples,'String',contents,'Value',1);

active_names = cellstr(get(H.Active_samples,'String'));

selected_rmv = active_selected{1,1};

for i = 1:length(selected_rmv)
	pdps_active(:,selected_rmv(1,i)) = 0;
end

pdps_active( :, all(~pdps_active,1) ) = [];

x = x_min:x_int:x_max;

axes(H.stacked_sink_samples);
nsamples = length(pdps_active(1,:));
if nsamples == 0
	cla(H.stacked_sink_samples,'reset');
elseif nsamples == 1
	plot(x,pdps_active, 'k','linewidth',1.5);
	axis([x_min x_max 0 max(max(pdps_active))])
else 
	for i = 1:nsamples
		dist_max(i,1) = max(pdps_active(:,i));
	end
	for i = 1:nsamples
		dist_max2(i,1) = sum(dist_max(1:i,1));
	end
	pdp_sum = sum(dist_max2);
	for i = 1:nsamples
		pdp_adj(:,i) = pdps_active(:,i) - dist_max2(i,1) + dist_max2(end,1);
	end
	for i = 1:nsamples
		plot(x,pdp_adj(:,i), 'k','linewidth',1.5);
	hold on
	end
	axis([x_min x_max 0 dist_max2(end,1)])
end

xlabel('Age (Ma)')
ylabel('Relative Probability')
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])

%Mean Kuiper
if nsamples==1;
    set(H.mean_V_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_V_sink(i,j)=max(cumsum(pdps_active(:,i))-...
                cumsum(pdps_active(:,j)))+max(cumsum(pdps_active(:,j))...
                -cumsum(pdps_active(:,i)));
        end
    end
    mean_V_sink_temp=[];
    for i=1:nsamples
        mean_V_sink_temp=vertcat(mean_V_sink_temp,diag(mean_V_sink,i));
    end
    mean_V_sink=round(mean(mean_V_sink_temp),3);
    range_V_sink=round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end

set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
guidata(hObject, H);


if isempty(active_names) == 1
	set(H.xaxis_min,'Enable','on');
	set(H.xaxis_max,'Enable','on');
	set(H.xaxis_int,'Enable','on');
	set(H.pdps,'Enable','on');
	set(H.kdes,'Enable','on');
elseif length(char(active_names)) == 0
	set(H.xaxis_min,'Enable','on');
	set(H.xaxis_max,'Enable','on');
	set(H.xaxis_int,'Enable','on');
	set(H.pdps,'Enable','on');
	set(H.kdes,'Enable','on');
end

end

%% PUSHBUTTON Merge and Activate
function Merge_Callback(hObject, eventdata, H)
opengl hardware
pdps_active = H.pdps_active;
text = H.text;
N = H.N;
data = H.data;
loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_selected = H.loaded_selected;

if isempty(loaded_names) == 1
	err_dlg=errordlg('There are no samples to merge!','Hang on...');
	waitfor(err_dlg);
elseif length(char(loaded_names)) == 0
	err_dlg=errordlg('There are no samples to merge!','Hang on...');
	waitfor(err_dlg);
else
	
	set(H.xaxis_min,'Enable','off');
	set(H.xaxis_max,'Enable','off');
	set(H.xaxis_int,'Enable','off');
	set(H.pdps,'Enable','off');
	set(H.kdes,'Enable','off');

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));


tmp = loaded_names(loaded_selected{1,1});
result = strcat(tmp(1:end-1,1), ',');
result(end+1,1) = tmp(end,1);
namecat = result{1,1};
for i = 2:length(result)
	namecat = [namecat, result{i,1}];
end

active_names = cellstr(get(H.Active_samples,'String'));
active_names(all(cellfun('isempty',active_names),2),:) = [];
active_names = [active_names; namecat];
set(H.Active_samples, 'String', active_names);

contents = get(H.Loaded_samples,'String');
if length(contents)<1; return; end %if already empty, do nothing
Index=get(H.Loaded_samples,'Value');
contents(Index)=[]; %remove the item
Value=Index-1;
if Value<1; Value=1;end %take care of exception
set(H.Loaded_samples,'String',contents,'Value',1);

comb = (strsplit(active_names{end,1},','))';

for i = 1:N
	for j = 1:length(comb)
		if strcmp(text(1,i*2-1), comb(j,1)) == 1
			comb_data_name{:,i*2-1} = comb(j,1);
			comb_data(:,i*2-1:i*2) = data(:,i*2-1:i*2);
		end
	end
end

comb_data_name = comb_data_name(~cellfun('isempty',comb_data_name));
comb_data(isnan(comb_data))=0;
comb_data( :, all(~comb_data,1) ) = [];

for i = 1:length(comb_data_name)
	comb_data_length(i,1) = length(nonzeros(comb_data(:,i*2-1)));
end

merged_data = [];

for i = 1:length(comb_data_name)
	merged_data(1+length(merged_data):length(merged_data)+comb_data_length(i,1),1:2) = comb_data(1:comb_data_length(i,1),i*2-1:i*2);
end

x = x_min:x_int:x_max;

if get(H.pdps,'Value') == 1
	m = nonzeros(merged_data(:,1));
	m = m(isfinite(m(:,1)),:);
	s = nonzeros(merged_data(:,2));
	s = s(isfinite(s(:,1)),:);
	f = zeros(length(m),length(x));
	for j = 1:length(m)
		f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*x_int);
	end
	pdps_merged(:,1) = ((sum(f, 1))/length(m)).';
end

if get(H.kdes,'Value') == 1
	rad_on_kernel=get(H.Dist_opts,'selectedobject');
	switch rad_on_kernel
		case H.kern_opt
			xA = transpose(x);
			dist_data = nonzeros(merged_data(:,1));
			[bandwidth,kdeA,xmesh1,cdf]=kde(dist_data,length(x),x_min,x_max);
			pdps_merged(:,1) = transpose(interp1(xmesh1, kdeA, xA));
		case H.kern_set
			dist_data = nonzeros(merged_data(:,1));
			kernel = str2num(get(H.kern,'String'));
			kernel_dist_data(1:length(nonzeros(merged_data(:,1))),1) = kernel;
			pdps_merged(:,1) = pdp5(dist_data(:,1),kernel_dist_data,x_min,x_max,x_int); 
	end
end


pdps_active = [pdps_active,pdps_merged];

axes(H.stacked_sink_samples);
nsamples = length(pdps_active(1,:));
if nsamples == 1
	plot(x,pdps_active, 'k','linewidth',1.5);
	axis([x_min x_max 0 max(max(pdps_active))])
else
	for i = 1:nsamples
		dist_max(i,1) = max(pdps_active(:,i));
	end
	for i = 1:nsamples
		dist_max2(i,1) = sum(dist_max(1:i,1));
	end
	pdp_sum = sum(dist_max2);
	for i = 1:nsamples
		pdp_adj(:,i) = pdps_active(:,i) - dist_max2(i,1) + dist_max2(end,1);
	end
	for i = 1:nsamples
		plot(x,pdp_adj(:,i), 'k','linewidth',1.5);
	hold on
	end
	axis([x_min x_max 0 dist_max2(end,1)])
end

xlabel('Age (Ma)')
ylabel('Relative Probability')
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])
%Mean Kuiper
if nsamples==1;
    set(H.mean_V_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_V_sink(i,j)=max(cumsum(pdps_active(:,i))-...
                cumsum(pdps_active(:,j)))+max(cumsum(pdps_active(:,j))...
                -cumsum(pdps_active(:,i)));
        end
    end
    mean_V_sink_temp=[];
    for i=1:nsamples
        mean_V_sink_temp=vertcat(mean_V_sink_temp,diag(mean_V_sink,i));
    end
    mean_V_sink=round(mean(mean_V_sink_temp),3);
    range_V_sink=round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end

set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
guidata(hObject, H);

end

%% PUSHBUTTON Activate All
function Activate_all_Callback(hObject, eventdata, H)
opengl hardware
name = H.name;
N = H.N;
text = H.text;
data = H.data;
loaded_names = cellstr(get(H.Loaded_samples,'String'));

if isempty(loaded_names) == 1
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
elseif length(char(loaded_names)) == 0
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
else
	
	set(H.xaxis_min,'Enable','off');
	set(H.xaxis_max,'Enable','off');
	set(H.xaxis_int,'Enable','off');
	set(H.pdps,'Enable','off');
	set(H.kdes,'Enable','off');

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

set(H.Active_samples, 'String', name);
set(H.Loaded_samples, 'String', []);
set(H.Loaded_samples, 'Value',1);

active_names = cellstr(get(H.Active_samples,'String'));
active_names(all(cellfun('isempty',active_names),2),:) = [];

for i = 1:N
	for j = 1:length(active_names)
		if strcmp(text(1,i*2-1), active_names(j,1)) == 1
			active_data_name{:,i*2-1} = active_names(j,1);
			active_data(:,i*2-1:i*2) = data(:,i*2-1:i*2);
		end
	end
end

active_data(isnan(active_data))=0;
active_data( :, all(~active_data,1) ) = [];
active_data_name = active_data_name(~cellfun('isempty',active_data_name));

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));

x = x_min:x_int:x_max;

if get(H.pdps,'Value') == 1
	for i = 1:length(active_data_name)
		m = nonzeros(active_data(:,i*2-1));
		m = m(isfinite(m(:,1)),:);
		s = nonzeros(active_data(:,i*2));
		s = s(isfinite(s(:,1)),:);
		f = zeros(length(m),length(x));
		for j = 1:length(m)
			f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*x_int);
		end
		pdps_active(:,i) = ((sum(f, 1))/length(m)).';
	end
end

if get(H.kdes,'Value') == 1
	rad_on_kernel=get(H.Dist_opts,'selectedobject');
	switch rad_on_kernel
		case H.kern_opt
			xA = transpose(x);
			for i = 1:length(active_data_name)
				dist_data = nonzeros(active_data(:,i*2-1));
				[bandwidth,kdeA,xmesh1,cdf]=kde(dist_data,length(x),x_min,x_max);
				pdps_active(:,i) = transpose(interp1(xmesh1, kdeA, xA));
				clear dist_data kdeA
			end
		case H.kern_set
			for i = 1:length(active_data_name)
				dist_data = nonzeros(active_data(:,i*2-1));
				kernel = str2num(get(H.kern,'String'));
				kernel_dist_data(1:length(nonzeros(active_data(:,i*2-1))),1) = kernel;
				pdps_active(:,i) = pdp5(dist_data(:,1),kernel_dist_data,x_min,x_max,x_int); 
			end
		end
end

axes(H.stacked_sink_samples);
nsamples = length(active_data_name);
if nsamples == 1
	plot(x,pdps_active, 'k','linewidth',1.5);
	axis([x_min x_max 0 max(max(pdps_active))])
else
	for i = 1:nsamples
		dist_max(i,1) = max(pdps_active(:,i));
	end
	for i = 1:nsamples
		dist_max2(i,1) = sum(dist_max(1:i,1));
	end
	pdp_sum = sum(dist_max2);
	for i = 1:nsamples
		pdp_adj(:,i) = pdps_active(:,i) - dist_max2(i,1) + dist_max2(end,1);
	end
	for i = 1:nsamples
		plot(x,pdp_adj(:,i), 'k','linewidth',1.5);
	hold on
	end
	axis([x_min x_max 0 dist_max2(end,1)])
end

xlabel('Age (Ma)')
ylabel('Relative Probability')
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])

H.pdps_active = pdps_active;

%Mean Kuiper
if nsamples==1;
    set(H.mean_V_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_V_sink(i,j)=max(cumsum(pdps_active(:,i))-...
                cumsum(pdps_active(:,j)))+max(cumsum(pdps_active(:,j))...
                -cumsum(pdps_active(:,i)));
        end
    end
    mean_V_sink_temp=[];
    for i=1:nsamples
        mean_V_sink_temp=vertcat(mean_V_sink_temp,diag(mean_V_sink,i));
    end
    mean_V_sink=round(mean(mean_V_sink_temp),3);
    range_V_sink=round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end

set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
guidata(hObject, H);

end

%% PUSHBUTTON Deactivate All.
function Deactivate_all_Callback(hObject, eventdata, H)

active_names = cellstr(get(H.Active_samples,'String'));

if isempty(active_names) == 1
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
elseif length(char(active_names)) == 0
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
else

	set(H.xaxis_min,'Enable','on');
	set(H.xaxis_max,'Enable','on');
	set(H.xaxis_int,'Enable','on');
	set(H.pdps,'Enable','on');
	set(H.kdes,'Enable','on');
	
cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

pdps_active = [];
name = H.name;

set(H.Active_samples, 'String', []);
set(H.Loaded_samples, 'String', name);
set(H.Active_samples, 'Value',1);

ylabel([])
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])

set(H.nsamples, 'String', 'N/A');
set(H.mean_R2_sink, 'String', 'N/A');
set(H.range_R2_sink, 'String', 'N/A');
set(H.mean_V_sink, 'String', 'N/A');
set(H.range_V_sink, 'String', 'N/A');
set(H.num_sources_opt,'String',20);

active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
guidata(hObject, H);

end

%% PUSHBUTTON Replot Results %%
function replot_res_Callback(hObject, eventdata, H)

optimized = H.optimized;

if 	optimized == 0
	err_dlg=errordlg('You need to determine the optimal number of sources or upload previous results.','Wait a sec...');
	waitfor(err_dlg);
elseif optimized == 1 && get(H.check_nmf_opt, 'Value') == 1
	xdata = H.xdata;
	ydata = H.ydata;
	SSR = H.SSR;
	modelpred = H.modelpred;

	figure;
	plot(xdata,ydata,'o',sort(xdata),modelpred,'r-');
	title('Final Residual versus Rank')
	xlabel('Rank')
	ylabel('Final Residual')

	figure;
	plot(xdata,SSR,'-o');
	title('SSR versus Rank')
	xlabel('Rank')
	ylabel('SSR1 + SSR2')
	
	optimized = 1;
	H.optimized = optimized;
	H.xdata = xdata;
	H.ydata = ydata;
	H.SSR = SSR;
	guidata(hObject, H);
end


%% PUSHBUTTON Export Plots %%
function export_plots_Callback(hObject, eventdata, H)

optimized = H.optimized;
nsources = H.nsources;
pdps_active = H.pdps_active;
source_PDP = H.source_PDP;

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));

x = x_min:x_int:x_max;


figure;
nsamples = length(pdps_active(1,:));
if nsamples == 1
	plot(x,pdps_active, 'k','linewidth',1.5);
	axis([x_min x_max 0 max(max(pdps_active))])
else
	for i = 1:nsamples
		dist_max(i,1) = max(pdps_active(:,i));
	end
	for i = 1:nsamples
		dist_max2(i,1) = sum(dist_max(1:i,1));
	end
	pdp_sum = sum(dist_max2);
	for i = 1:nsamples
		pdp_adj(:,i) = pdps_active(:,i) - dist_max2(i,1) + dist_max2(end,1);
	end
	for i = 1:nsamples
		plot(x,pdp_adj(:,i), 'k','linewidth',1.5);
	hold on
	end
	axis([x_min x_max 0 dist_max2(end,1)])
end

xlabel('Age (Ma)')
ylabel('Relative Probability')
%ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])
title('Sink Samples')

clear dist_max2

figure;
if nsources == 1
	plot(x,source_PDP, 'k','linewidth',1.5);
else
	for i = 1:nsources
		dist_max(i,1) = max(source_PDP(:,i));
	end
	for i = 1:nsources
		dist_max2(i,1) = sum(dist_max(1:i,1));
	end
	pdp_sum = sum(dist_max2);
	for i = 1:nsources
		pdp_adj(:,i) = source_PDP(:,i) - dist_max2(i,1) + dist_max2(end,1);
	end
	for i = 1:nsources
		plot(x,pdp_adj(:,i), 'k','linewidth',1.5);
	hold on
	end
	axis([x_min x_max 0 dist_max2(end,1)])
end

xlabel('Age (Ma)')

ylabel('Relative Probability')
%ytickformat('%.0f')
set(H.stacked_reconstructed_sources,'ytick',[]);
axis([x_min x_max 0 dist_max2(end,1)])
title('Reconstructed Sources')


if optimized == 1 

	xdata = H.xdata;
	ydata = H.ydata;
	SSR = H.SSR;
	modelpred = H.modelpred;

	figure;
	plot(xdata,ydata,'o',sort(xdata),modelpred,'r-');
	title('Final Residual versus Rank')
	xlabel('Rank')
	ylabel('Final Residual')

	figure;
	plot(xdata,SSR,'-o');
	title('SSR versus Rank')
	xlabel('Rank')
	ylabel('SSR1 + SSR2')
	
	optimized = 1;

	H.xdata = xdata;
	H.ydata = ydata;
	H.SSR = SSR;
	
	
end

H.optimized = optimized;
guidata(hObject, H);

%% PUHBUTTON Example Input %%
function pushbutton20_Callback(hObject, eventdata, H)

Example_Data_Set


















































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








































%% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


%% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)



%% PUSHBUTTON Export NMF Results
%% --- Executes on button press in Export_NMF_Results.
function Export_NMF_Results_Callback(hObject, eventdata, H)
global iter_num;
nsources = H.nsources;
export = H.export;

if get(H.check_nmf, 'Value') == 1
	BaseName=' Sources';
	tst3{1,1}=[num2str(nsources),BaseName];
	[file,path] = uiputfile('*.xls','Save file');
	%warning('off','MATLAB:xlswrite:AddSheet');
	if ismac == 1
		xlwrite([path file], export, char(tst3(1,1)));
	end
	if ispc == 1
		xlswrite([path file], export, char(tst3(1,1)));
	end
	%RemoveSheet123([path file]);
	%warning('off');
end

if get(H.check_nmf_opt, 'Value') == 1	
	tst3{1,iter_num-1} = [];
	BaseName=' Sources';
	for i=2:iter_num
		tst3{1,i-1}=[num2str(i),BaseName];
	end
	[file,path] = uiputfile('*.xls','Save file');
	warning('off','MATLAB:xlswrite:AddSheet');
	if ismac == 1
		for i=1:iter_num-1
			xlwrite([path file], export(:,:,i), char(tst3(1,i)));
		end
	end
	if ispc == 1
		for i=1:iter_num-1
			xlswrite([path file], export(:,:,i), char(tst3(1,i)));
		end
	end	
	
	
	%RemoveSheet123([path file]);
	%warning('off');
end







	
	
	
	

function iterations_Callback(hObject, eventdata, H)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterations as text
%        str2double(get(hObject,'String')) returns contents of iterations as a double


%% --- Executes during object creation, after setting all properties.
function iterations_CreateFcn(hObject, eventdata, H)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



































%% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


%% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)





%% --- Executes during object creation, after setting all properties.
function Active_samples_CreateFcn(hObject, eventdata, H)
% hObject    handle to Active_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









function num_sources_Callback(hObject, eventdata, H)
% hObject    handle to num_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_sources as text
%        str2double(get(hObject,'String')) returns contents of num_sources as a double


%% --- Executes during object creation, after setting all properties.
function num_sources_CreateFcn(hObject, eventdata, H)
% hObject    handle to num_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Final_Residual_Callback(hObject, eventdata, H)
% hObject    handle to Final_Residual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Final_Residual as text
%        str2double(get(hObject,'String')) returns contents of Final_Residual as a double


%% --- Executes during object creation, after setting all properties.
function Final_Residual_CreateFcn(hObject, eventdata, H)
% hObject    handle to Final_Residual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tof_Callback(hObject, eventdata, H)
% hObject    handle to tof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tof as text
%        str2double(get(hObject,'String')) returns contents of tof as a double


%% --- Executes during object creation, after setting all properties.
function tof_CreateFcn(hObject, eventdata, H)
% hObject    handle to tof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, H)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, H)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)












function num_sources_opt_Callback(hObject, eventdata, H)
% hObject    handle to num_sources_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_sources_opt as text
%        str2double(get(hObject,'String')) returns contents of num_sources_opt as a double


% --- Executes during object creation, after setting all properties.
function num_sources_opt_CreateFcn(hObject, eventdata, H)
% hObject    handle to num_sources_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)






function xaxis_min_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_min as text
%        str2double(get(hObject,'String')) returns contents of xaxis_min as a double


% --- Executes during object creation, after setting all properties.
function xaxis_min_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xaxis_max_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_max as text
%        str2double(get(hObject,'String')) returns contents of xaxis_max as a double


% --- Executes during object creation, after setting all properties.
function xaxis_max_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xaxis_int_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_int as text
%        str2double(get(hObject,'String')) returns contents of xaxis_int as a double


% --- Executes during object creation, after setting all properties.
function xaxis_int_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kern_opt.
function kern_opt_Callback(hObject, eventdata, H)
% hObject    handle to kern_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kern_opt


% --- Executes on button press in kern_set.
function kern_set_Callback(hObject, eventdata, H)
% hObject    handle to kern_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kern_set



function kern_Callback(hObject, eventdata, H)
% hObject    handle to kern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kern as text
%        str2double(get(hObject,'String')) returns contents of kern as a double


% --- Executes during object creation, after setting all properties.
function kern_CreateFcn(hObject, eventdata, H)
% hObject    handle to kern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function opt_num_result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opt_num_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called





% --- Executes on button press in check_nmf_opt.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to check_nmf_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_nmf_opt


% --- Executes on button press in check_nmf.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to check_nmf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_nmf


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% PUSHBUTTON Cancel %%
function cancel_Callback(hObject, eventdata, handles)

global cancel
cancel = 1;
