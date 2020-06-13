function varargout = ProductionPlanning(varargin)
% PRODUCTIONPLANNING MATLAB code for ProductionPlanning.fig
%      PRODUCTIONPLANNING, by itself, creates a new PRODUCTIONPLANNING or raises the existing
%      singleton*.
%
%      H = PRODUCTIONPLANNING returns the handle to a new PRODUCTIONPLANNING or the handle to
%      the existing singleton*.
%
%      PRODUCTIONPLANNING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRODUCTIONPLANNING.M with the given input arguments.
%
%      PRODUCTIONPLANNING('Property','Value',...) creates a new PRODUCTIONPLANNING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProductionPlanning_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProductionPlanning_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProductionPlanning

% Last Modified by GUIDE v2.5 25-Aug-2019 15:07:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ProductionPlanning_OpeningFcn, ...
    'gui_OutputFcn',  @ProductionPlanning_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before ProductionPlanning is made visible.
function ProductionPlanning_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProductionPlanning (see VARARGIN)

% Choose default command line output for ProductionPlanning
handles.output = hObject;
set(handles.radiobutton2,'Value',1);
set(handles.initguess_check,'Value',0);
set(handles.init_guess,'Visible','off');
set(handles.termination_state,'Visible','off');
set(handles.table1,'Visible','off')
set(handles.table2,'Visible','off')
set(handles.save_dir,'enable','on');
set(handles.save_dir,'Visible','off');
set(handles.save_browse,'Visible','off');
set(handles.ppp,'enable','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ProductionPlanning wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ProductionPlanning_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function budget_Callback(hObject, eventdata, handles)
% hObject    handle to budget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of budget as text
%        str2double(get(hObject,'String')) returns contents of budget as a double


% --- Executes during object creation, after setting all properties.
function budget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to budget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rawmat_avail_Callback(hObject, eventdata, handles)
% hObject    handle to rawmat_avail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rawmat_avail as text
%        str2double(get(hObject,'String')) returns contents of rawmat_avail as a double

% --- Executes during object creation, after setting all properties.
function rawmat_avail_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rawmat_avail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nproducts_Callback(hObject, eventdata, handles)
% hObject    handle to nproducts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nproducts as text
%        str2double(get(hObject,'String')) returns contents of nproducts as a double


% --- Executes during object creation, after setting all properties.
function nproducts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nproducts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_time_Callback(hObject, eventdata, handles)
% hObject    handle to max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_time as text
%        str2double(get(hObject,'String')) returns contents of max_time as a double


% --- Executes during object creation, after setting all properties.
function max_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_file.
function browse_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ExPath
[FileName1,FilePath1]=uigetfile('.xls');
Path = [FilePath1 FileName1];
if any(Path)~=0
    ExPath = Path;
    set(handles.edit9,'String',ExPath);
else
    ExPath = get(handles.edit9,'String');
end


function init_guess_Callback(hObject, eventdata, handles)
% hObject    handle to init_guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of init_guess as text
%        str2double(get(hObject,'String')) returns contents of init_guess as a double
global initial_guess;
initial_guess = str2num(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function init_guess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to init_guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in start_push.
function start_push_Callback(hObject, eventdata, handles)
% hObject    handle to start_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ExPath handles stopit bp_val
stopit=0;
set(handles.nprocesses,'enable','off');
set(handles.nproducts,'enable','off');
set(handles.budget,'enable','off');
set(handles.rawmat_avail,'enable','off');
set(handles.edit9,'enable','off');
set(handles.max_time,'enable','off');
set(handles.bp_check,'enable','off');
set(handles.bp_check,'enable','off');
set(handles.termination_state,'Visible','on');
set(handles.termination_state,'String','intlinprog is running');
set(handles.reset_push,'Enable','off');
set(handles.stop_push,'Enable','on');
set(handles.save_dir,'enable','off');
set(handles.ppp,'enable','off');
set(handles.init_guess,'enable','off');
set(handles.init_guess_browse,'enable','off');

y = xlsread(ExPath,'data');
[~,~,yy] =  xlsread(ExPath,'data');

product=y(:,2);

l=y(:,3);
m=y(:,4);
h=y(:,5);
il=y(:,6);
im=y(:,7);
ih=y(:,8);
cl=y(:,9);
cm=y(:,10);
ch=y(:,11);
saleprice=y(:,12);

R = str2num(handles.rawmat_avail.String);
R = R(:);
B = str2num(handles.budget.String);


nr = length(R);     % number of raw materials
np = length(il);    % number of processes

consumeraw = y(:,13:12+nr);


nprocesses = str2num(handles.nprocesses.String);

if bp_val == 0
    z = [];
    NbyP  = 0;
else
    z = xlsread(ExPath,'byproducts');
    [~,~,zz] =  xlsread(ExPath,'byproducts');
    Sjk = [];
    byProduct = z(:,1);% No of by products by each process
    maxP = max(byProduct);
    
    SelP = z(:,2*maxP+2:3*maxP+1);
    Alpha = z(:,2:maxP+1);
    Beta = z(:,maxP+2:2*maxP+1);
    NbyP  = 0;ubBP = [];
    for i = 1:nprocesses
        NbyP = NbyP + byProduct(i);
        Sjk = [Sjk; SelP(i,1:byProduct(i))'];
        ubBP = [ubBP repmat(h(i),1,byProduct(i))];
    end
end





% x = [X L M H Y Z]; Arrangement of decision variables



lb = zeros(6*np,1); % lower bound zeros
ub = [h;ones(5*np,1)]; % upper bounds for continuous and binary variables

if bp_val == 1
    lb = [lb;zeros(NbyP,1)];
    ub = [ub;ubBP'];
end
    
    
intcon = 4*np+1:6*np; % location of the binary variables

ZM = zeros(np,np); % A matrix of zeros that is often required in the model

%Coeefficients of the objective function
Profit = [saleprice; -cl; -cm; -ch; -zeros(2*np,1)];
if bp_val==1
    Profit = [Profit;Sjk];
end
A1 = []; A2 = []; A3 = []; A4 = []; A5 = []; b = [];
A1 = [ZM eye(np) ZM ZM -eye(np) ZM]; % Equation 6
A2 = [ZM ZM ZM eye(np) eye(np) ZM]; % Equation 7
A3 = [eye(np) ZM ZM ZM ZM -10^6*eye(np)]; % Equation 9
A4 = [zeros(1,np) il' im' ih' zeros(1,2*np)]; % Equation 13
A5 = [consumeraw' zeros(nr,5*np)]; % raw materials Equation 11

% Combining all the matrices of the inequality constraints

uni = handles.radiobutton5.Value;
if uni == 1 || handles.radiobutton6.Value
    nt = max(product);  % number of products
    UniMat = zeros(nt,np);
    
    for p = 1:np
        UniMat(product(p),p) = 1;
    end
    A6 = [zeros(nt,5*np),UniMat]; % unique process constraint
    A = [A1;A2;A3;A4;A5;A6];
    if handles.radiobutton6.Value
        b = [zeros(np,1);ones(np,1);zeros(np,1);B;R;str2num(handles.ppp.String)';];
    else
        b = [zeros(np,1);ones(np,1);zeros(np,1);B;R;ones(nt,1)];
    end
else
    A = [A1;A2;A3;A4;A5];
    b = [zeros(np,1);ones(np,1);zeros(np,1);B;R];
end

if bp_val==1
    A = [A zeros(size(A,1),NbyP)];
end
    

% Combining all the matrices of the equality constraints
A1 = [eye(np) -l.*eye(np) -m.*eye(np) -h.*eye(np) ZM ZM ];
A2 = [ZM eye(np) eye(np) eye(np) ZM -eye(np)];
if bp_val==1
    A1 = [A1,zeros(np,NbyP)];
    A2 = [A2,zeros(np,NbyP)];
end
Aeq = [A1;A2];
beq = zeros(2*np,1);

%% by-product
if bp_val==1
    Aeq2 = [];
    QzV = eye(NbyP);
    count = 0;
    for j = 1:np
        nKj = byProduct(j);
        XMAT = zeros(nKj,np);
        XMAT(:,j) = 1;
        bZM = zeros(nKj,4*np);
        alpha = Alpha(j,1:byProduct(j));
        beta = Beta(j,1:byProduct(j));
        %     QzV(count+1:count+nKj,count+1:count+nKj) = eye(nKj);
        Aeq2 = [Aeq2;  -alpha'.*XMAT   bZM  -beta'.*XMAT QzV(count+1:count+nKj,:)];
        count = count+nKj;
    end
    beq2 = zeros(NbyP,1);
    
    Aeq = [Aeq;Aeq2];
    beq = [beq;beq2];
end

% Solving the MILP problem
global initial_guess
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxTime',str2num(handles.max_time.String),'OutputFcn',@stopmyiter,'PlotFcn',@optimplotmilp);
[x,fval,exitflag,output] = intlinprog(-Profit,intcon,A,b,Aeq,beq,lb,ub,initial_guess,options);
saveas(gcf,'Convergence_curve.png')
set(handles.reset_push,'Enable','on');
set(handles.stop_push,'Enable','off');

switch exitflag
    case 2
        exitval = 'Status: intlinprog stopped prematurely. Integer feasible point found.';
    case 1
        exitval = 'Status: intlinprog converged to the solution x.';
    case 0
        exitval = 'Status: intlinprog stopped prematurely. No integer feasible point found.';
    case -1
        exitval = 'Status: intlinprog stopped by an output function or plot function.';
    case -2
        exitval = 'Status: No feasible point found.';
    case -3
        exitval = 'Status: Root LP problem is unbounded.';
end
set(handles.termination_state,'Visible','on');
set(handles.termination_state,'String',exitval);
if ~stopit

    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALL FUNCTION FILE
[R_utilized, B_utilized, product_used, prodprofit,product_produced, process_used, procprofit,amt_produced]  = ppanalysis(x,l, m, h, cl, cm, ch, A4, A5, saleprice,product,NbyP);
   
   
axes(handles.axes2);
pie(prodprofit,product_used);
title('Profit breakdown per product');

axes(handles.axes3);
pie(procprofit,process_used);
title('Profit breakdown per process');

figure1 = figure;
pie(prodprofit,product_used);
title('Profit breakdown per product');
saveas(figure1,'PPT.png') 


figure2 = figure;
pie(procprofit,process_used);
title('Profit breakdown per process');
saveas(figure2,'PPP.png')


handles.table1.RowName = process_used;
handles.table1.ColumnName = 'Amt Produced';
handles.table1.Data = amt_produced;


handles.table2.RowName = product_used;
handles.table2.ColumnName = 'Amt Produced';
product_produced = product_produced(product_produced~=0);
handles.table2.Data = product_produced';
set(handles.table1,'Visible','on')
set(handles.table2,'Visible','on')


set(handles.budget_util,'String',num2str(B_utilized));
set(handles.rm_util,'String',num2str([round(R_utilized,1)']));
set(handles.net_profit,'String',num2str(abs(fval)));
%% Byprod plot
if bp_val 
    byProdAmt = x(6*np+1:end);
    maxByProd = max(byProduct);
    byProdMat = NaN(np,maxByProd);
    byProdMat(1,1:byProduct(1)) = byProdAmt(1:byProduct(1));
    count = byProduct(1);
    for p = 2:np
        byProdMat(p,1:byProduct(p)) = ...
            byProdAmt(count+1:count+byProduct(p));
        count = count+byProduct(p);
    end
    byProdMat(byProdMat==0) = NaN;
    hold off
    figure
    plot(1:np,byProdMat,'*','MarkerSize',8);
    xlabel('Process number')
    ylabel('Amount of by-products produced')
    legend('q_1','q_2','q_3','q_4','q_5','Location','Best')
    saveas(gcf,'BP Produced.png')
end

global DirPath
check_save = handles.save_check.Value;
if check_save ==1
    clk = clock;
    x_save = x(1:np);
    ind = x_save>0;
    activeP = find (x_save>0);
    x_save = x_save(activeP);
    l_save = x(np+1:2*np);
    l_save = l_save(ind);
    m_save = x(2*np+1:3*np);
    m_save = m_save(ind);
    h_save = x(3*np+1:4*np);
    h_save = h_save(ind);
%     y_save = x(4*np+1:5*np);
%     z_save = x(5*np+1:6*np);
    vec_save = [x_save l_save m_save h_save];
    cell_save = num2cell(vec_save);
    if bp_val
        head_save = {'X' 'L' 'M' 'H','q_1','q_2','q_3','q_4','q_5'};
        producedByProd = num2cell(byProdMat(activeP,:));
        save_val1 = [head_save;cell_save producedByProd];
    else
        head_save = {'X' 'L' 'M' 'H'};
        save_val1 = [head_save;cell_save];
    end
    
    name = [DirPath '/' num2str(clk(4)) num2str(clk(5)) num2str(ceil(clk(6)))];
    head_fval = {'Fval'};
    save_fval = [head_fval,num2cell(abs(fval))];
    xlswrite(name,save_fval,'results','B1:C1');
    xlswrite(name,['Active process';process_used'],'results',strcat('A2:A',num2str(length(x_save)+2)));
    if bp_val
        xlswrite(name,save_val1,'results',strcat('B2:J',num2str(length(x_save)+2)));
    else
        xlswrite(name,save_val1,'results',strcat('B2:E',num2str(length(x_save)+2)));
    end
    xlswrite(name,yy,'main_data');
    
    if bp_val==1
        xlswrite(name,zz,'bp_data');
    end
end      
end

% --- Executes on button press in reset_push.
function reset_push_Callback(hObject, eventdata, handles)
% hObject    handle to reset_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stopit bp_val
set(handles.init_guess,'Visible','off');
set(handles.init_guess_browse,'Visible','off');
set(handles.init_guess,'enable','on');
set(handles.init_guess_browse,'enable','on');
set(handles.nprocesses,'enable','on','String','54');
set(handles.nproducts,'enable','on','String','24');
set(handles.budget,'enable','on','String','2000');
set(handles.rawmat_avail,'enable','on','String','[1000 1000]');
set(handles.edit9,'enable','on','String','');
set(handles.max_time,'enable','on','String','7200');
set(handles.init_guess,'enable','off');
set(handles.ppp,'enable','off');
%set(handles.ppp,'enable','on','String','');
set(handles.termination_state,'Visible','off');
set(handles.initguess_check,'Value',0);
set(handles.init_guess,'Visible','off');
set(handles.save_check,'enable','on','Value',0);
set(handles.save_dir,'enable','on','String','');
set(handles.save_dir,'Visible','off');
axes(handles.axes2);
y = gca; cla reset;
y.Title.String = '';
handles.table2.RowName = '';
handles.table2.ColumnName = '';
handles.table2.Data = [];
handles.table1.RowName = '';
handles.table1.ColumnName = '';
handles.table1.Data = [];
set(handles.table1,'Visible','off')
set(handles.table2,'Visible','off')
set(handles.budget_util,'enable','on','String','');
set(handles.rm_util,'enable','on','String','');
set(handles.net_profit,'enable','on','String','');
axes(handles.axes3); cla reset;
y = gca;
y.Title.String = '';

set(handles.bp_check,'enable','on','Value',0);
stopit = 0;
bp_val=0;


% --- Executes on button press in stop_push.
function stop_push_Callback(hObject, eventdata, handles)
% hObject    handle to stop_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stopit
stopit=1;
set(handles.reset_push,'Enable','on');

% --- Executes on button press in help_push.
function help_push_Callback(hObject, eventdata, handles)
% hObject    handle to help_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open([pwd '\Help_doc.pdf'])

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_check.
function save_check_Callback(hObject, eventdata, handles)
% hObject    handle to save_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_check
if handles.save_check.Value
    set(handles.save_dir,'Visible','on');
    set(handles.save_browse,'Visible','on');
else
    set(handles.save_dir,'Visible','off')
    set(handles.save_browse,'Visible','off');
end

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton2.Value
    set(handles.radiobutton5,'Value',0)
    set(handles.radiobutton6,'Value',0)
    set(handles.ppp,'enable','off');
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton2



function nprocesses_Callback(hObject, eventdata, handles)
% hObject    handle to nprocesses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nprocesses as text
%        str2double(get(hObject,'String')) returns contents of nprocesses as a double


% --- Executes during object creation, after setting all properties.
function nprocesses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nprocesses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function termination_state_Callback(hObject, eventdata, handles)
% hObject    handle to termination_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of termination_state as text
%        str2double(get(hObject,'String')) returns contents of termination_state as a double


% --- Executes during object creation, after setting all properties.
function termination_state_CreateFcn(hObject, eventdata, handles)
% hObject    handle to termination_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in formulation_push.
function formulation_push_Callback(hObject, eventdata, handles)
% hObject    handle to formulation_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open([pwd '\formulation_gui_1.pdf'])

% --- Executes when entered data in editable cell(s) in table1.
function table1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to budget_util (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of budget_util as text
%        str2double(get(hObject,'String')) returns contents of budget_util as a double


% --- Executes during object creation, after setting all properties.
function budget_util_CreateFcn(hObject, eventdata, handles)
% hObject    handle to budget_util (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to rm_util (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rm_util as text
%        str2double(get(hObject,'String')) returns contents of rm_util as a double


% --- Executes during object creation, after setting all properties.
function rm_util_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rm_util (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton5.Value
    set(handles.radiobutton2,'Value',0)
    set(handles.radiobutton6,'Value',0)
    set(handles.ppp,'enable','off');
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton6.Value
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton5,'Value',0);
    set(handles.ppp,'visible','on');
    set(handles.ppp_browse,'visible','on');
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton6


function ppp_Callback(hObject, eventdata, handles)
% hObject    handle to ppp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ppp as text
%        str2double(get(hObject,'String')) returns contents of ppp as a double


% --- Executes during object creation, after setting all properties.
function ppp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in initguess_check.
function initguess_check_Callback(hObject, eventdata, handles)
% hObject    handle to initguess_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of initguess_check
ig = handles.initguess_check.Value;
global initial_guess;
if ig ==1
    set(handles.init_guess,'Visible','on');
    set(handles.init_guess_browse,'Visible','on');
    set(handles.init_guess,'enable','on');
    set(handles.init_guess_browse,'enable','on');
    initial_guess = str2num(handles.init_guess.String);
else
    set(handles.init_guess,'visible','off');
    set(handles.init_guess_browse,'visible','off');
    set(handles.init_guess,'enable','off');
    set(handles.init_guess_browse,'enable','off');
    initial_guess = [];
end



function save_dir_Callback(hObject, eventdata, handles)
% hObject    handle to save_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_dir as text
%        str2double(get(hObject,'String')) returns contents of save_dir as a double


% --- Executes during object creation, after setting all properties.
function save_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_browse.
function save_browse_Callback(hObject, eventdata, handles)
% hObject    handle to save_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DirPath
path=uigetdir();
if any(path)~=0
    DirPath = path;
    set(handles.save_dir,'String',DirPath);
else
    DirPath = pwd;
end


% --- Executes on button press in bp_check.
function bp_check_Callback(hObject, eventdata, handles)
% hObject    handle to bp_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bp_check
global bp_val
bp_val = handles.bp_check.Value;


% --- Executes on button press in ppp_browse.
function ppp_browse_Callback(hObject, eventdata, handles)
% hObject    handle to ppp_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pppPath
[FileName1,FilePath1]=uigetfile('.xls');
Path = [FilePath1 FileName1];
if any(Path)~=0
    pppPath = Path;
    ppp_data = xlsread(pppPath);
    set(handles.ppp,'String',num2str(ppp_data'));
end

% --- Executes on button press in init_guess_browse.
function init_guess_browse_Callback(hObject, eventdata, handles)
% hObject    handle to init_guess_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global initgPath
[FileName1,FilePath1]=uigetfile('.xls');
Path = [FilePath1 FileName1];
if any(Path)~=0
    initgPath = Path;
    init_data = xlsread(initgPath);   
    set(handles.init_guess,'String',num2str(init_data'));
%    ExPath = get(handles.edit9,'String');
end