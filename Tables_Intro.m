clear
clc

% Creating, storing, reading, and manipulating tables which may contain
% arbitrary data types

% Get arrays of the same length - the variable name will be column headers
% and the value will be the column

% Start producing initial table
n = 9; % number of rows for our table
% Init numerical stuff
zerosbase = zeros(n,1);
H = zerosbase;
h = zerosbase;

% Init cell array for ARBITRARY DATA TYPES (can be vectors!)
% To access the contents of a cell, enclose indices in curly braces e.g.
% cellbase{1} = contents_of_first_entry
cellbase = cell(n,1); % like zeros but for cell
Avec = cellbase; % Avec = A vector

for i = 1:n
    H(i) = i;
    h(i) = i/10;

    % Init A vec for demonstration purposes of storing a vector
    A = zeros(1,n);
    for j = 1:n
        A(j) = i*sin((2*pi*(j-1))/n);
    end
    Avec{i} = A;

end

% Print each
disp('Input lists for our table:')
H
h
Avec
disp('First entry of A:')
Avec{1}
disp('Notice the above is a vector, but Avec itself shows as a bunch of 1xn matrices rather than the values')

% Tables are then defined using the following syntax, and produces a table
% with n rows (not including headers), n is the length of each vector, and
% the headers are the variable names respectively

disp('The table is as follows:')
T = table(H,h,Avec)
% To specify the column names manually can use
% T = table(H,h,Avec,'VariableNames',{'Column1name','Column2name','Column3name'})

% Can access columns using attribute notation - dots
% Here is accessing the column Avec, and getting the first entry
disp('Can access individual columns:')
disp(T.Avec{1})


% Will now define a new similar table, and demonstrate combining them
m = 2*n;

% Init numerical stuff
zerosbase = zeros(m,1);
H = zerosbase;
h = zerosbase;

% Init cell array for arbitrary data types (can be vectors!)
% To access the contents of a cell, enclose indices in curly braces e.g.
% cellbase{1} = contents_of_first_entry
cellbase = cell(m,1);
Avec = cellbase;

for i = 1:m
    H(i) = (i+max(1,n-3));
    h(i) = (i+max(1,n-3))/10;

    % Init A vec for demonstration purposes of storing a vector
    A = zeros(1,n);
    for j = 1:n
        A(j) = i*sin((2*pi*(j-1))/n);
    end
    Avec{i} = A;

end

% Tables are then defined using the following syntax, and produces a table
% with n rows (not including headers), n is the length of each vector, and
% the headers are the variable names respectively
disp('Second demo table with some new values:')
T2 = table(H,h,Avec)

% Get the rows in T2 which are not in T1 
% WARNING: prone to floating point errors!
% Taddnaive = setdiff(T2,T) % Does not work due to Avec being a weird data
% type

% Solution: 
% Define two new tables which have only the parameters for our
% table, in this case H and h which lie in column positions 1 and 2 of our
% table.
disp('Get first two columns of T (and T2) for setdifference:')
key_columns = [1,2] % e.g [2,4] represents 2nd and 4th column
Tparams = T(:,key_columns) % Like matrix notation
T2params = T2(:,key_columns);

% Can get the new parameters in T2 with setdiff, and find the row numbers
% for these new parameters (as a vector for the row numbers)
disp('New parameters and their T2 row positions:')
[Tnewparams,T2rowpositions] = setdiff(T2params,Tparams)

% Using similar notation to getting the columns we want, we can get the
% rows we want from T2 using our vector above
disp('Table including new parameters and their values from T2:')
T2new = T2(T2rowpositions,:)

% Combine T with T2new by directly defining a new table using block matrix
% notation
disp('Table with all the found parameters and their values')
T = [T; T2new]

% Write to a .mat file - this allows us to preserve things like the data
% types (ensures the special columns are read and written properly)
tablename = "tables_intro_table.mat";
save(tablename,"T",'-v7.3')
disp(['Written T to ',tablename])

% Overwrite T for demo purposes
T = 0;

% Read from the .mat file, test retrieving Avec (the problem child)
disp(['Read the variable T from ',tablename])
load(tablename,"T")
T
disp('First Avec: ')
T.Avec{1}

% Also demo writing the table to an xml file
tablenamexml = "tables_intro_table.xml";
writetable(T,tablenamexml)
disp(['Written T to ',tablenamexml])


% Demo combining a third table to the read table T
m = 2*m;


% Init numerical stuff
zerosbase = zeros(m,1);
H = zerosbase;
h = zerosbase;

% Init cell array for arbitrary data types (can be vectors!)
% To access the contents of a cell, enclose indices in curly braces e.g.
% cellbase{1} = contents_of_first_entry
cellbase = cell(m,1);
Avec = cellbase;

for i = 1:m
    H(i) = (i+max(1,n-3));
    h(i) = (i+max(1,n-3))/10;

    % Init A vec for demonstration purposes of storing a vector
    A = zeros(1,n);
    for j = 1:n
        A(j) = i*sin((2*pi*(j-1))/n);
    end
    Avec{i} = A;

end

% Tables are then defined using the following syntax, and produces a table
% with n rows (not including headers), n is the length of each vector, and
% the headers are the variable names respectively
disp('Third demo table with some new values:')
T2 = table(H,h,Avec)

% Get the rows in T2 which are not in T1 
% Define two new tables which have only the parameters for our
% table, in this case H and h which lie in column positions 1 and 2 of our
% table.
key_columns = [1,2];
Tparams = T(:,key_columns);
T2params = T2(:,key_columns);

% Can get the new parameters in T2 with setdiff, and find the row numbers
% for these new parameters (as a vector for the row numbers)
[Tnewparams,T2rowpositions] = setdiff(T2params,Tparams);

% Using similar notation to getting the columns we want, we can get the
% rows we want from T2 using our vector above
T2new = T2(T2rowpositions,:);

% Combine T with T2new by directly defining a new table using block matrix
% notation
disp('Table with all the found parameters and their values')
T = [T; T2new]

% Overwrite the tables in our .mat (and xml for good measure)
tablename = "tables_intro_table.mat";
save(tablename,"T",'-v7.3')
disp(['Written T to ',tablename])

tablenamexml = "tables_intro_table.xml";
writetable(T,tablenamexml)
disp(['Written T to ',tablenamexml])

