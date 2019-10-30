function [ level_set ] = BFS_level_set( Adj, s )
% Find level sets based on Breadth First Search 
%
% @Xiaozhe Hu, Department of Mathematics, Tufts University
% 01/03/2018
%
% Input:       
%       s:      starting vertices
%       Adj:     adjancy matrix
%
% Output: 
%       level_set:  level set      
%   

%---------------------
% initilize
%---------------------
n = size(Adj,1);        % get number of vertices

Q = s;                % initialize the queue and push s into it

level_set = zeros(n,1); % initilize level set
level_set(s) = 1;

% label = -ones(n,1);     % initialize the label
%                         % -1: unvisited
%                         %  0: visited
%                         %  1: explored

%---------------------
% main loop
%---------------------
while ~isempty(Q)
   
    v = Q(1);           % pop the first vetex in Q
    Q(1) = [];
        
    [j,~,~] = find(Adj(:,v));  % find neighbors of vetex v
    
    for w = j'           % loop over neighbors
        
        if level_set(w) == 0   % if it is unvisited
            
            level_set(w) = level_set(v)+1;   % mark it as visited
            Q = [Q, w];     % append it to queue
            
        end
        
    end
    
end

end

