function [ C, label, E ] = mDFS( Adj, s, label)
% Modified Depth First Search (use stack)
%
% @Xiaozhe Hu, Department of Mathematics, Tufts University
% 01/29/2016
%
% Input:       
%       s:      starting vertex
%       Adj:      adjancy matrix
%   label:      lable of vertices (-1: unvisited; 0: visited; >0:label for 2DFS) 
%
% Output: 
%       C:      connected component including vetex s
%       E:      egeds on the DFS tree
%   label:      updated label of vertices
%


%---------------------
% initilize
%---------------------
n = size(Adj,1);         % get number of vertices

C = [s];                 % connected component
E = [];                  % edges on DFS tree

j = 1;                   % label for 2DFS                  

S = [s];                 % initialize stack
label(s) = 0;            % label the starting vertex as visited

flag = true;             % flag of whether we can go further

while ~isempty(S)
   
    while flag
        
        v = S(1);                      % get the first vertex in stack
        
        [Nv,~,~] = find(Adj(:,v));       % find its neighbors
        
        idx = find(label(Nv) == -1);   % check any adjacent vertex is unvisited
        
        if ~isempty(idx)               % if there is such vertex
            
            w = Nv(idx(1));            % choose one of such vertex
            label(w) = 0;              % label it us visited
            C = [C,w];                 % put it in connected component
            E = [E, [v;w]];            % put the edges 
            S = [w, S];                % push onto stack
        
        else
            
            flag = false;              % otherwise, get out of the loop
            
        end
        
    end
    
    label(S(1)) = j;                   % before backtrack, label the vertex
    j = j+1;
    
    S(1) = [];                         % pop stack (backtrack)
    flag = true;                            
        
end