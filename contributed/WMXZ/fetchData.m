% Fetch seismic data from the internet
% WMXZ
% uses GISMO
% please edit the following functions to meet your requirements
% startup: needs path to GISMO library
% stationList: 
%     SL: provide a list of stations to be shown
%              statuon name, channel, network, location, source name
%     Server: source name, irisFetch flag, dtatSource or web url
%     filter object: 
%
function fetchData
    % must be first function, so check setup and initialization
    if ~exist('datasource','file')
        startup;
    end
    initFetch();
end

%%%%%>>>>>>>> begin user specific area (to be edited)
function startup()
    addpath('C:\Users\Walter\Documents\Github\GISMO') % path to GISMO
    startup_GISMO
end

function [SL,Server,s]=stationList()
% station list to be used
    SL={...{'RDF14', 'SHZ', 'AM', '00', 'AM'};
        ...{'GRA1',  'BHZ', 'GR', '*', 'BGR'};
        ...{'GRB5',  'BHZ', 'GR', '*', 'BGR'};
        ...{'GRC1',  'BHZ', 'GR', '*', 'BGR'};
        {'RA78B', 'SHZ', 'AM', '00', 'L2'};
        ...{'MSSA',  'BHZ', 'IV', '*', 'INGV'};
%        {'EQUI',  'HHZ', 'GU', '*', 'INGV'}
        };
    %
    % define data sources or web services for other than L1/2, AM and IRIS-DMC
    % L1/L2, AM and DMC are accessed via waveform method (2nd element is 0)
    % other stations are accessed via irisFetch (2nd element is 1)
    Server = {{'L1',  0, datasource('winston', '192.168.0.8', 16032)}; % local raspberryshakedata
              {'L2',  0, datasource('winston', '192.168.0.102', 16032)}; % local raspberryshakedata
              {'AM',  0, datasource('winston', 'raspberryshakedata.com', 16032)}; %RS community
              {'DMC', 0, datasource('irisdmcws')}; % data_source is IRIS DMC webservices
              {'BGR', 1, 'eida.bgr.de'};
              {'INGV',1, 'webservices.ingv.it'}};
    %
    % oject for common filter
    % parameters: fft size, overlap, LP corner frequency, min-max dB level
    s=spectralobject(1024, round(0.9*1024), 5, 40+[0 60]);
end
%%%%%<<<<<<<< end of user specific area
%
%------------------------------------------------------------------------
% no need to edit after this line (but yuo are free to do it)
function initFetch()
global f duration
    f = figure('Visible','off');
    set(f,'position',[40   60   960   900]);    % change dimension if needed

    duration = 3600;                            % change default duration if needed

    f.Name = ' Fetch Seismic Data (WMXZ)';      % change figure name if desired
    initGui();
    f.Visible='on';
    colormap(f,'jet');                          % change colomap if needed
end

function initGui()
global h3 startText duration

    if isempty(startText)
        startText = datestr(now(), 'dd-mm-yyyy HH:MM:SS');
    end

    clf
    uicontrol('Style','pushbutton',...
                 'String','now','Position',[20,870,70,25],...
                 'Callback',@next_Callback,'fontsize',12);
    uicontrol('Style','pushbutton',...
                 'String','load','Position',[120,870,70,25],...
                 'Callback',@load_Callback,'fontsize',12);

    uicontrol('Style','text','String','Start:',...
            'HorizontalAlignment','right',...
           'Position',[230,870-2,70,25],'fontsize',12);
    h3  = uicontrol('Style','edit','String',startText,...
           'Position',[300,870,200,25],...
           'Callback',@edit_Callback,'fontsize',12);

    uicontrol('Style','text','String','Duration:',...
            'HorizontalAlignment','right',...
           'Position',[500,870-2,100,25],'fontsize',12);
    uicontrol('Style','edit','String',duration,...
           'Position',[600,870,100,25],...
           'Callback',@duration_Callback,'fontsize',12);

end

function next_Callback(source,eventdata) 
   getData(0);
end

function load_Callback(source,eventdata) 
    getData(1);
end

function edit_Callback(source,eventdata) 
global startText
    startText=get(source,'String');
end

function duration_Callback(source,eventdata) 
global duration
    duration=str2num(get(source,'String'));
end

function set_Timestamp(txt)
global h3
        set(h3,'string',txt);
end

function w = getWave(scnl,startTime,endTime,web)
    tr=irisFetch.Traces(get(scnl,'network'), ...
                        get(scnl,'station'),...
                        get(scnl,'location'),...
                        get(scnl,'channel'),startTime,endTime,web);
    if length(tr)>1, tr=tr(1); end
    w = waveform(scnl, tr.sampleRate, tr.startTime, tr.data);
end

function getData(type)
global f startText duration

    oldpointer = get(f, 'pointer'); 
    oldpointer = 'arrow'; % uncommented if testing and prog crashed 
    set(f, 'pointer', 'watch') 
    drawnow;

    [SL,Server,s]=stationList();
    
    NSL=length(SL);

    scnl = repmat(scnlobject('','','',''),NSL,1);
    for ii=1:NSL
        scnl(ii) = scnlobject(SL{ii}{1},SL{ii}{2},SL{ii}{3},SL{ii}{4});
    end

    if type== 0
        t_end=(now()-1/24)-0/(24*60);
        t_start=t_end-duration/(24*3600);
        startTime = datestr(t_start,'yyyy-mm-dd HH:MM:SS');
        endTime   = datestr(t_end,  'yyyy-mm-dd HH:MM:SS');
        startText = datestr(t_start,'dd-mm-yyyy HH:MM:SS');
        set_Timestamp(startText);
    else
        t_piv=datenum(startText,'dd-mm-yyyy HH:MM:SS');
        startTime = datestr( t_piv, 'yyyy-mm-dd HH:MM:SS');
        endTime   = datestr( t_piv + duration/(24*3600), 'yyyy-mm-dd HH:MM:SS');
    end

    fprintf('Selected times are:\n%s\n%s\n', startTime,endTime);

    w=repmat(waveform,NSL,1);
    for ii=1:NSL
        source = SL{ii}{5};
        fprintf('Requesting Data for %s ...\n',SL{ii}{1});
        %
        ds=[];
%        web='service.iris.edu';
        web=[];
        for jj=1:length(Server)
            if strcmp(source,Server{jj}{1}) % find always last entry
                if Server{jj}{2}==0
                    ds=Server{jj}{3};
                elseif Server{jj}{2}==1
                    web = Server{jj}{3};  
                end
            end
        end
        %
        if ~isempty(ds)
            w(ii) = waveform(ds, scnl(ii), startTime, endTime);
        elseif ~isempty(web)
            w(ii) = getWave(scnl(ii),startTime,endTime,['http://', web]);
        else
            fprintf('no datasource or web service defined for %s\n',source);
            w(ii) = waveform(scnl(ii), 0, startTime, []); % empty trace
        end
    end
    
    fprintf('Flter data ...\n');

    fobj = filterobject('b', [0.2 5], 2);
    wx=[];
    for ii=1:length(w)
        if ~isempty(get(w(ii),'data')) 
            w(ii) = fillgaps(w(ii),'meanall');
            wx = [wx,filtfilt(fobj, w(ii))]; 
        end
    end

    fprintf('Display data ...\n');
    %
    initGui();
    plotData(wx,s);
    f.Visible = 'on';
    set(f, 'pointer', oldpointer)

end

function plotData(w,s)
global f
    nw=length(w);

    hh=0.85/nw;
    ax=[];
    for ii=1:nw
        yo=0.1+(ii-1)*hh;
        ax(2*ii-1)=axes('position',[0.1 yo,0.8,hh*0.75],'fontsize',12);
        yo=0.1+(ii-1)*hh +hh*0.75;
        ax(2*ii)=axes('position',[0.1 yo,0.8,hh*0.25],'fontsize',12);
    end

    nfft = round(get(s,'nfft'));
    overlap = floor(get(s, 'over'));
    dBlims = get(s, 'dBlims');
    fmax = get(s, 'freqmax');

    maxT=0;
    minT=inf;
    for ii=1:nw
        fsamp = get(w(ii), 'freq');
        data = get(w(ii), 'data');
        To = get(w(ii),'start');
        minT = min(minT, To);
        maxT = max(maxT,To+length(data)/fsamp/(24*3600));
    end
    minT=mod(minT,1);
    maxT=mod(maxT,1);
    if maxT<minT, minT=minT-1; end
    xl=[minT,maxT];
    
    for ii=1:nw
        jj=1+nw-ii; % plot first station to the top
        fsamp = get(w(jj), 'freq');
        data = get(w(jj), 'data');
        To = get(w(jj),'start');
        
        if fsamp>40, nfftx=nfft; else nfftx=nfft/2; end
        tt=mod(To + (0:length(data)-1)'/fsamp/(24*3600),1);
        if tt(1)>tt(end), tt(tt>tt(end))=tt(tt>tt(end))-1; end
            
        overlap = round(0.9*nfftx);

        [S,F,Ts] = spectrogram(data, nfftx, overlap, nfftx, fsamp);
        Y = 20*log10(abs(S)+eps);

        T=mod(To + Ts/(24*3600),1); 
        if T(1)>T(end), T(T>T(end))=T(T>T(end))-1; end 
        index = find(F <= fmax);
        F=F(index);
        Y=Y(index,:);

        i1=2*ii-1;
        i2=2*ii;
        
        plot(ax(i2),tt,data)
        imagesc(T,F,Y,'parent',ax(i1),dBlims);
        %
        set(ax(i1),'fontsize',12)
        
        axis(ax(i1),'xy');

        ylabel(ax(i1), sprintf('%s\n%s\n%s',...
            get(w(jj), 'network'),...
            get(w(jj), 'station'),...
            get(w(jj), 'channel')), 'FontSize', 12);
        xlabel(ax(i1),'')
        set(ax(i2),'ytick',[]);
        set(ax(i1),'xlim',xl);
        set(ax(i2),'xlim',xl);
        %
        if ii>1
            set(ax(i1),'xtick',get(ax(1),'xtick'));
            set(ax(i2),'xtick',get(ax(1),'xtick'));
            set(ax(i1),'xticklabel',[]);
        else
            datetick(ax(1),'x','HH:MM','keeplimits')
            set(ax(i2),'xtick',get(ax(1),'xtick'));
        end
        set(ax(i2),'xticklabel',[]);
    end
    linkaxes(ax,'x')
    h=zoom(f);
    setAxesZoomMotion(h,ax,'horizontal');
    fprintf('%f %f\n',mean(data),std(data));
end