list2maple := proc(ls :: list(nonnegint));
    local i;
    return add( seq( ls[i] * x^(i - 1), i = 1 .. nops(ls) ) );
end proc;


verifyFile := proc(in_name :: string, to_check :: posint, out_name :: string, $)
    local infile,outfile,i,a,b,c,p,res;
    
    p := 2^16 + 1;
    infile := fopen(in_name, READ);
    outfile := fopen(out_name, WRITE);

    for i to to_check do
        a := list2maple(parse(readline(infile)));
        b := list2maple(parse(readline(infile)));
        c := list2maple(parse(readline(infile)));
        res := Expand(a*b) mod p;
        if res <> c then
            fprintf(outfile, "%a\n", [ seq( coeff(res, x, i), i=1..degree(res) ) ]);
        fi;
    od;

    fclose(infile);
    fclose(outfile);
end proc;
