p := 2^16 + 1;


benchmark := proc(filename :: string, $)
    local file,size,total_time,a,b,t1,t2,res;

    file := fopen("benchmark", WRITE);

    for size from 9 to 15 do
        a := Randpoly(2^size, x) mod p;
        b := Randpoly(2^size, x) mod p;
        t1 := Time();
        res := Expand(a*b) mod p;
        t2 := Time();
        total_time := Value(t2) - Value(t1);
        fprintf(file, "%d\n", total_time);
    od;

    fclose(file);
end proc;
