function [project_workspace_path] = Initialization_Env

[~, machine_name] = system('hostname');
machine_name = cellstr(machine_name);

% switch machine_name{1}
%     case 'DESKTOP-5R0S9NF'  % Old computer before replace SDD
%         project_workspace_path = 'D:/CRCNS_Project/';
%     case 'DESKTOP-RBIV19A'  % New computer after replace SDD
%         project_workspace_path = 'D:/CRCNS_Project/';
%     case 'VPN-128-237-147-85.LIBRARY.VPN.CMU.EDU'
%         project_workspace_path = '/Users/chenyu/Workspace/CRCNS Project/';
%     case 'Chens-MacBook-Pro.local'
%         project_workspace_path = '/Users/chenyu/Workspace/CRCNS Project/';
%     case 'Chens-MBP.fios-router.home'
%         project_workspace_path = '/Users/chenyu/Workspace/CRCNS Project/';
%     case [char(10) 'psych-o.hpc1.cs.cmu.edu']
%         project_workspace_path = '/Users/chenyu/Workspace/CRCNS Project/';
%     case 'psych-o.hpc1.cs.cmu.edu'
%         project_workspace_path = '/Users/chenyu/Workspace/CRCNS Project/';
%     case 'Chens-MBP.wv.cc.cmu.edu'
%         project_workspace_path = '/Users/chenyu/Workspace/CRCNS Project/';
%     otherwise
%         disp('This is a new environment, please figure out workspace path first!')
%         return

%% judge using new method

if ispc
    project_workspace_path = 'D:/CRCNS_Project/';
elseif ismac
    project_workspace_path = '/Users/chenyu/Workspace/CRCNS Project/';
else
    disp('This is a new environment, please figure out workspace path first!')
    return    

end
