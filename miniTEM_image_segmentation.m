%Copyright (c) 2019 Andrea Behanova, Ali Abdollazadeh
% A.I. Virtanen Inbstitute for Molecular Sciences, University of Eastern
% Finland, Finland

%version 1.0  13.12.2019

%Permission is hereby granted, free of charge, to any person obtaining a copy 
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.

% This software uses Bio-Formats 5.9.2 package, which can be found in https://www.openmicroscopy.org/bio-formats/downloads/

function varargout = miniTEM_image_segmentation(varargin)
% MINITEM_IMAGE_SEGMENTATION MATLAB code for miniTEM_image_segmentation.fig
%      MINITEM_IMAGE_SEGMENTATION, by itself, creates a new MINITEM_IMAGE_SEGMENTATION or raises the existing
%      singleton*.
%
%      H = MINITEM_IMAGE_SEGMENTATION returns the handle to a new MINITEM_IMAGE_SEGMENTATION or the handle to
%      the existing singleton*.
%
%      MINITEM_IMAGE_SEGMENTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MINITEM_IMAGE_SEGMENTATION.M with the given input arguments.
%
%      MINITEM_IMAGE_SEGMENTATION('Property','Value',...) creates a new MINITEM_IMAGE_SEGMENTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before miniTEM_image_segmentation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to miniTEM_image_segmentation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help miniTEM_image_segmentation

% Last Modified by GUIDE v2.5 13-Dec-2019 15:14:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @miniTEM_image_segmentation_OpeningFcn, ...
                   'gui_OutputFcn',  @miniTEM_image_segmentation_OutputFcn, ...
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


% --- Executes just before miniTEM_image_segmentation is made visible.
function miniTEM_image_segmentation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to miniTEM_image_segmentation (see VARARGIN)

currentFolder = pwd;
addpath([currentFolder '\bfmatlab'])

handles.data = [];
handles.indx = [];
handles.file = [];
handles.virus =[];

% Choose default command line output for miniTEM_image_segmentation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes miniTEM_image_segmentation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = miniTEM_image_segmentation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function open_file_Callback(hObject, eventdata, handles)
% hObject    handle to open_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Files, pathname]=uigetfile('*.ome.tif', 'Chose files to load:','MultiSelect','on');

if size(Files,2) == 1
    return
end
    

if iscell(Files)
    L = size(Files,2);
    f = waitbar(1/L,'Loading files');
else
    L = 1;
    filename = Files;
end

for i = 1:L
    if iscell(Files)
        filename = Files{1,i};
    end
    data = bfopen([pathname filename]);
    s1 = data{1,1}{1,1};
    %save to the memory
    handles.data = [handles.data;{s1}];

    %pop up
    txt = get(handles.file_popup,'String');
    txt = [txt; {filename}];
    set(handles.file_popup,'String',txt);


    % Adding file to file list
    txt = get(handles.Files_listbox,'String');
    txt = [txt; {filename}];
    set(handles.Files_listbox,'String',txt);

    files = get(handles.Files_listbox,'String');
    handles.file = files(2:end,1);
    
    if iscell(Files)
        waitbar(i/L,f)
    end
    
end

if iscell(Files)
    close(f)
end


guidata(hObject, handles);


% --- Executes on selection change in Files_listbox.
function Files_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Files_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Files_listbox


% --- Executes during object creation, after setting all properties.
function Files_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in file_popup.
function file_popup_Callback(hObject, eventdata, handles)
% hObject    handle to file_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns file_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from file_popup
contents = cellstr(get(hObject,'String'));
display_choice = contents(get(hObject, 'Value'));
handles.dispChoice = display_choice;
handles.virus = [];

if strcmp(handles.dispChoice, 'Select file...')
    set(get(gca,'children'),'Visible','off')
    handles.dispChoice = [];
elseif isempty(handles.dispChoice)
    return
else
    indx = strcmp(handles.Files_listbox.String, handles.dispChoice);
    indx = indx(2:end);
    handles.indx = indx;
    data = handles.data{indx};
    
    imshow(data)
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function file_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in segmentation_button.
function segmentation_button_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADENOVIRUSES SEGMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(handles.dispChoice)
    return
end

indx = handles.indx;
data = handles.data{indx};

%% Histogram adjustment
Imad = imadjust(data);

%% Median filtering
J = medfilt2(Imad,[15 15]);


%% Selection of low intensity areas
th = round(max(J(:))/1.5);
J(J>th)=max(J(:));
J(J<=th)=0;

J=~logical(J);
J = imdilate(J,true(15));
J = imfill(J,'holes');

J = bwareaopen(J,2000);

Image = data;
Image(J==0) = 0;
Imad(J==0) = 0;

%% Hough Transform
[centers, radii, metric] = imfindcircles(Imad,[16 50],'EdgeThreshold',0.0001,'ObjectPolarity','bright','Sensitivity',0.87);
centers = round(centers);

Imad = imadjust(data);

%% Removing candidates od adenoviruses which do not have dar area around

if length(centers)>70 %before70  %if lot of FP change it to lower
    for i = 1:length(centers)
        
        x1 = centers(i,2)-100;
        x2 = centers(i,2)+100;
        y1 = centers(i,1)-100;
        y2 = centers(i,1)+100;

        if x1<1
            x1=1;
        end
        
        if y1<1
            y1=1;
        end
        
        if x2>size(Imad,2)
            x2=size(Imad,2);
        end
        
        if y2>size(Imad,1)
            y2=size(Imad,1);
        end

        candidate = Imad(x1:x2,y1:y2);

        h = imhist(candidate);
        [int,peak] = max(h);
        if peak >50
            centers(i,1)=0;
            centers(i,2)=0;
            radii(i)=0;
        end
    end
end

%% Visualization

numbers = find(centers(:,1));
centers = centers(numbers,:);
radii = nonzeros(radii);

mask = zeros(size(Imad));

for i = 1:size(centers,1)
    circ = drawcircle('Center',centers(i,:),'Radius',radii(i));
    BW = createMask(circ);
    delete(circ)
    mask(BW)=i;
end



handles.virus.centers = centers;
handles.virus.radii = radii;
handles.virus.mask = mask;

guidata(hObject, handles);

% h = fspecial('sobel');
% blurred = imfilter(J,h,'replicate'); 
% imshow(blurred);

% BW1 = edge(J,'Sobel');
% imshow(BW1)
Disp(hObject, eventdata, handles)

function Disp(hObject, eventdata, handles)
indx = handles.indx;
data = handles.data{indx};

Imad = imadjust(data);
imshow(Imad)

centers = handles.virus.centers;
radii = handles.virus.radii;

viscircles(centers, radii,'EdgeColor','r');

% figure
% imshow(Imad);

% --- Executes on button press in add_label_button.
function add_label_button_Callback(hObject, eventdata, handles)
% hObject    handle to add_label_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.virus)
    return
end

indx = handles.indx;
data = handles.data{indx};

centers = handles.virus.centers;
radii = handles.virus.radii;
mask = handles.virus.mask;

h = imellipse;
pos = round(getPosition(h));
BW = createMask(h);
delete(h)
mask(BW)=max(mask(:))+1;
center = round([pos(1)+pos(3)/2, pos(2)+pos(4)/2]);
radius = pos(3)/2;
centers = [centers;center];
radii = [radii; radius];

handles.virus.centers = centers;
handles.virus.radii = radii;
handles.virus.mask = mask;

guidata(hObject, handles);
Disp(hObject, eventdata, handles)


% --- Executes on button press in remove_label.
function remove_label_Callback(hObject, eventdata, handles)
% hObject    handle to remove_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.virus)
    return
end

centers = handles.virus.centers;
radii = handles.virus.radii;
mask = handles.virus.mask;

indx = handles.indx;
data = handles.data{indx};


[y,x] = getpts; y = round(y(:,1)); x = round(x(:,1));
index = [];
for i = 1:length(x)
    if (x(i)>0) && (x(i)<size(data,2)) && (y(i)>0) && (y(i)<size(data,1))
        linearInd = sub2ind([size(data,2),size(data,1)], x(i), y(i));
        index(i) = nonzeros(mask(linearInd));
        if index(i)==max(mask(:))
            centers(index(i),:) = [];
            radii(index(i)) = [];
        else
            centers(index(i),:) = 0;
            radii(index(i)) = 0;
        end
        mask(mask==index(i))=0;
          
    else
        return
    end
end


handles.virus.centers = centers;
handles.virus.radii = radii;
handles.virus.mask = mask;


guidata(hObject, handles);
Disp(hObject, eventdata, handles)

% --- Executes on button press in save_label_button.
function save_label_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_label_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.virus)
    return
end

centers = handles.virus.centers;
radii = handles.virus.radii;
mask = handles.virus.mask;

indx = handles.indx;
name = handles.file{indx};
name = [name(1:end-8), '_target.mat'];
save(name,'mask', '-v7.3')


% --- Executes on button press in segmentation_debris.
function segmentation_debris_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation_debris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEBRIS SEGMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(handles.dispChoice)
    return
end

f = waitbar(0.1,'Debris Segmenttaion in process...');

indx = handles.indx;
data = handles.data{indx};

%% Histogram adjustment

Imad = imadjust(data);

%% Maximally stable extremal regions features detection

regions = detectMSERFeatures(Imad,'RegionAreaRange',[50 1500], 'ThresholdDelta',1);


mask = zeros(size(Imad));
for i = 1:length(regions.PixelList)
    ind = sub2ind(size(Imad),regions.PixelList{i,1}(:,2),regions.PixelList{i,1}(:,1));
    mask(ind)=1;
end

waitbar(0.2,f)


%% Filling holes in the segmentation
mask = imfill(mask,'holes');

%% Division of too close segments
m1 = imerode(mask,true(3));
m2 = bwlabel(m1,4);
mask = imdilate(m2,true(3));

waitbar(0.3,f)

%% Loading the mask of adenoviruses segmentation
[file,path] = uigetfile('*.mat','Load adenovirus mask');
adenovir = importdata(fullfile(path,file));

%% Removing adenoviruses
mask(logical(adenovir))=0;


%% Segmentation od rods
stats = regionprops(mask,'Area','MajorAxisLength','MinorAxisLength');

maj = [stats.MajorAxisLength];
min = [stats.MinorAxisLength];

non_Cidx = (maj>60);
non_Cidx2 = (min<30);

non_C = mask;

for i = 1:length(non_Cidx)
    if non_Cidx(i)==0 || non_Cidx2(i)==0 
        non_C(non_C==i)=0;
    end
end
rods = non_C;

mask(logical(rods)) = 0;

waitbar(0.5,f)

%% Division of the debris into smaller segments 

stats = regionprops(mask,'Area');
areas = cat(1,stats.Area);
standev = std(areas);

big_Aidx = (areas>standev);

big_A = mask;

for i = 1:length(big_Aidx)
    if big_Aidx(i)==0
        big_A(big_A==i)=0;
    end
end

waitbar(0.6,f)

I = imgaussfilt(Imad,2);

reg_max_big_A = imregionalmax(I);
reg_max_big_A(big_A == 0)=0;

D = bwdist(reg_max_big_A);
DL = watershed(D);

mask2 = logical(mask);
mask2(DL==0)=0;
mask2 = bwareaopen(mask2, 50);


m1 = imerode(mask2,true(3));
m2 = bwlabel(m1,4);
big_Debris = imdilate(m2,true(3));

%% Possible clustering testing

% stats = regionprops(big_Debris,'Area','MajorAxisLength','MinorAxisLength','Eccentricity','Solidity','Perimeter');
% areas = cat(1,stats.Area);
% major = cat(1,stats.MajorAxisLength);
% minor = cat(1,stats.MinorAxisLength);
% ecc = cat(1,stats.Eccentricity);
% sol = cat(1,stats.Solidity);
% allPerimeters = [stats.Perimeter];
% allAreas = [stats.Area];
% circ = (4*allAreas*pi)./(allPerimeters.^2);
% 
% vectors = [normalize(areas)'; normalize(major)'; normalize(minor)'; normalize(ecc)'; normalize(sol)'; normalize(circ)];
% 
% for i = 1:4
%     for j = 2:5
%         for k = 3:6
%             figure
%             scatter3(vectors(i,:),vectors(j,:),vectors(k,:))
%             switch i
%                 case 1
%                     x = 'Areas';
%                 case 2
%                     x = 'Major Axis Length';
%                 case 3
%                     x = 'Minor Axis Length';
%                 case 4
%                     x = 'Eccentricity';
%             end
%             switch j
%                 case 5
%                     y = 'Solidity';
%                 case 2
%                     y = 'Major Axis Length';
%                 case 3
%                     y = 'Minor Axis Length';
%                 case 4
%                     y = 'Eccentricity';
%             end
%             switch k
%                 case 5
%                     z = 'Solidity';
%                 case 6
%                     z = 'Circularity';
%                 case 3
%                     z = 'Minor Axis Length';
%                 case 4
%                     z = 'Eccentricity';
%             end
%             
%             
%             
%             xlabel(x)
%             ylabel(y)
%             zlabel(z)
%         end
%     end
% end

%% Division of debris into 2 groups: small and big debris based on a area threshold

stats = regionprops(big_Debris,'Area');
areas = cat(1,stats.Area);

big_Debris_idx = (areas>=200);

waitbar(0.7,f)

for i = 1:length(big_Debris_idx)
    if big_Debris_idx(i)==1
        big_Debris(big_Debris==i)=2;
    else
       big_Debris(big_Debris==i)=1; 
    end
end

waitbar(0.9,f)

mask2 = big_Debris;

%% Finalizing the mask

%add rods
mask2(logical(rods))=3;
mask2(logical(adenovir))=5;

%edges
J = imerode(mask2,true(6));
mask2 = mask2-J;

waitbar(1,f)

%% Visualization
axes(handles.Display)
imshow(Imad); hold on;
Lrgb = label2rgb(mask2,'lines', 'k');
himage = imshow(Lrgb); himage.AlphaData = 0.5;

close(f)
