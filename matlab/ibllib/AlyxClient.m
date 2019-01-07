classdef AlyxClient
    properties (SetAccess=private)
        base_url = ''
        user = ''
        timeout = 30
        weboptions = []
    end
    
    properties (SetAccess=private, Hidden=true)
        password = ''
        token = ''
        headers = ''
    end
    
    methods
        function self = AlyxClient(varargin)
            % ac = AlyxClient()
            % ac = AlyxClient('user','test','password','pass','base_url',...
            %                 'https://test.alyx.internationalbrainlab.org');
            prefs = self.getparameters(length(varargin)<6);
            % Handle input arguments, input arguments always overwrite preferences
            p = inputParser;
            addParameter(p,'user', prefs.user, @isstr)
            addParameter(p,'password', prefs.password, @isstr)
            addParameter(p,'base_url', prefs.base_url, @isstr)
            parse(p,varargin{:});
            for fn = fieldnames(p.Results)'; eval(['self.' fn{1} '= p.Results.' (fn{1}) ';']); end
            if isempty(self.password), self.password = prefs.password; end
            if isempty(self.base_url), self.base_url = prefs.base_url; end
            if isempty(self.user)    , self.user     = prefs.user;     end
            % setup weboptions for REST queries
            self.weboptions = weboptions(...
                'MediaType','application/json',...
                'Timeout',self.timeout, ...
                'CertificateFilename',''); %R2016b does not handle certificates well
            self = self.authenticate();
        end
    end
    
    methods (Access=private)
        function self = authenticate(self)
            % REST query to authenticate against Alyx and get an access token
            try
                rep = self.post('/auth-token', struct('username', self.user, 'password', self.password));
            catch ME
                error(['Connection issue while connecting to Alyx. Check your credentials !' char(10) ME.message])
            end
            self.token = rep.token;
            self.weboptions.HeaderFields = { 'Authorization', ['Token ' self.token]};
        end
        
        function prefs = getparameters(self, prompt)
            % Get parameters from preferences
            if nargin==0, prompt = true; end
            prefs = getpref('Alyx');
            if isempty(prefs)
                if ~prompt
                    prefs = struct('base_url','','user','','password','');
                    return
                end
                self.setup(prompt);
                prefs = getpref('Alyx');
            end
        end
        
        function url = format_url(self,url)
            if isempty(strfind(url, self.base_url))
                if ~startsWith(url, '/'), url = ['/' url]; end
                url = [self.base_url  url];
            end
        end
    end
    
    methods (Access = public)        
         function rep = get(self,url)
             % rep = get(url)
             % rep = ac.get('/sessions/86e27228-8708-48d8-96ed-9aa61ab951db')
             % rep = ac.get('https://test.alyx.internationalbrainlab.org/sessions/86e27228-8708-48d8-96ed-9aa61ab951db')
            rep = webread(self.format_url(url), self.weboptions);
            rep = flatten(rep);
         end
         
         function session_info = get_session(self, session_url)
             % session_info = ac.get_session('86e27228-8708-48d8-96ed-9aa61ab951db')
             % session_info = ac.get_session('https://test.alyx.internationalbrainlab.org/sessions/86e27228-8708-48d8-96ed-9aa61ab951db') 
            if isempty(strfind(session_url, self.base_url))
                session_url = [self.base_url '/sessions/' session_url];
            end
            % query the specific endpoint as the details as a slightly different output
            is = find(session_url=='/',1,'last');
            session_url = [session_url(1:is-1) '?id=' session_url(is+1:end) ];
            session_info = self.get(session_url);
         end
         
         function rep = post(self, url , request_struct)
             % rep = ac.post(url, request_struct)
            rep = webwrite(self.format_url(url),  jsonencode(request_struct),...
                setfield(self.weboptions, 'RequestMethod', 'post') );
         end
         
         function rep = put(self, url, request_struct)
             % rep = ac.put(url, request_struct)
            rep = webwrite(self.format_url(url),  jsonencode(request_struct),...
                setfield(self.weboptions, 'RequestMethod', 'put') );
         end
         
         function rep = delete(self, url)
             % rep = ac.delete(url)
            rep = webread(self.format_url(url),...
                setfield(self.weboptions, 'RequestMethod', 'delete'));
         end
%          function create_session(self, session_structure)
%              % self.create_session(session_structure)
%             %  session =  struct with fields: 
%             %        subject: 'clns0730'
%             %     procedures: {'Behavior training/tasks'}
%             %      narrative: 'auto-generated session'
%             %     start_time: '2018-07-30T12:00:00'
%             %           type: 'Base'
%             %         number: '1'
%             %          users: {'olivier'}
%          end
    end
    
        
    methods (Static)
        function setup()
            % AlyxClient.setup()
            % Prompts the user for base_url, user and password and stores for subsequent uses.
            prefs = getpref('Alyx');
            if isempty(prefs)
                prefs = struct('base_url','https://test.alyx.internationalbrainlab.org',...
                               'user','test_user',...
                               'password','');
            end
            if ~prompt, return, end
            % prompt for address
            base_url = input(['Alyx full URL: (example: https://test.alyx.internationalbrainlab.org), (current: ' prefs.base_url ') '], 's');
            if ~isempty(base_url)
                prefs.base_url = base_url;
            end
            % prompts for user
            user = input(['Alyx username (example:test_user), (current: ' prefs.user ')'], 's'); 
            if ~isempty(user)
                prefs.user = user;
            end
            % prompts for password
            % prefs.password
            password =  gui_password('DialogTitle', 'Enter password for Alyx instance');
            if ~isempty(password)
                prefs.password = password;
            end
            % assign properties
            for ff = fields(prefs)'
                setpref('Alyx', ff{1},  prefs.(ff{1}));
            end
        end
    end
end
