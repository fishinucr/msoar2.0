const char* defs_zmm[] = {
"",
"\\set{reqver}{2004-168}",
"\\if{\\cmp{gq}{\\__version__}{\\reqver}}{}{",
"\\write{stderr}{txt}{Requires zoem version >= \\reqver\\|}\\exit}",
"",
"\\${txt}{\\special{",
"   {-1}{ }",
"   {-2}{\\\\N}",
"   {-3}{-}",
"}}",
"\\${html}{\\special{",
"   {60}  {&lt;}      \\: less than sign",
"   {62}  {&gt;}      \\: greater than sign",
"   {38}  {&amp;}     \\: ampersand",
"",
"   {-1} {&nbsp;}     \\: the zoem escape \\~",
"   {-2} {<br>}       \\: the zoem escape \\|",
"   {-3} {-}          \\: the zoem escape \\-",
"}}",
"",
"\\${txt}{",
"   \\set{lref#2}{\\2}",
"   \\set{enref#2}{\\2}",
"   \\set{iref#2}{\\2}",
"   \\set{aref#2}{\\2}",
"   \\set{par}{\\@{\\P}}",
"   \\set{nl}{\\@{\\N}}",
"}",
"\\${html}{",
"   \\set{lref#2}{\\<a href=\"\\1\">{\\2}}",
"   \\set{iref#2}{\\<a href=\"#\\1\">{\\2}}",
"   \\set{aref#2}{\\<a href=\"\\1\">{\\2}}",
"   \\set{enref#2}{\\<a name=\"\\1\">{\\2}}",
"   \\set{par}{\\@{\\P}}",
"   \\set{nl}{\\@{\\N}}",
"}",
"\\${html}{",
"\\set{fmt_header}{",
"\\<html>",
"\\<head>",
"\\<style type=\"text/css\">",
"body { font-family: courier, mono; background: white }",
"a:link { text-decoration: none; color: #1111dd; }",
"a:visited { text-decoration: none; color: #11dd11; }",
"a:active { text-decoration: none; color: #dd1111; }",
"div.c { margin-left:6%; whitespace:pre; }",
"div.e { margin-left:6%; whitespace:pre; }",
"div.d { margin-left:12%; whitespace:pre; }",
"\\</style>",
"\\</head>",
"\\<body>",
"\\<pre>",
"}}",
"\\${txt}{\\set{fmt_header}{}}",
"",
"\\${html}{",
"   \\set{fmt_rule}{\\par\\<hr noshade size=1!>\\par\\${txt}{\\par}}",
"   \\set{fmt_footer}{\\</pre>\\</body>\\</html>}",
"}",
"\\${txt}{",
"   \\set{fmt_rule}{\\par\\format{%72~{-}{}.}{{}}\\par}",
"   \\set{fmt_footer}{}",
"}",
"",
"\\vanish{",
"",
"\\set{index_top}{\\enref{top}{}}",
"",
"   root info",
"   clustering{sz}{self}{cvg}{cvg}",
"",
"\\set{clustering#4}{\\enref{top}{}\\",
"Clustering has \\1 clusters\\nl",
"General measures: cohesion \\2, cov \\3, maxcov \\4}",
"",
"",
"   index clustering caption",
"\\set{cixcaption}{\\par\\enref{cl}{Clusters}\\nl}",
"\\set{cixend}{}",
"",
"   index element caption",
"\\set{eixcaption}{\\par\\enref{el}{Nodes}\\nl}",
"\\set{eixend}{}",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   Index links to node part, cluster part, and other index.",
"",
"                           \\formatted{",
"                           \\${txt}{",
"",
"\\set{ref_index_top#1}{}",
"",
"                           }\\${html}{",
"",
"\\set{ref_index_top#1}{",
"   \\enref{top}{}",
"   \\iref{cl}{Clusters}\\~",
"   \\iref{el}{Nodes}\\~",
"   \\iref{\\1.html}{\\1.html}",
"}",
"",
"                           }}",
"",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   cluster index entry",
"   \\\\cix{cl}{fname}{sz}{coh}{cov}",
"",
"                           \\formatted{",
"                           \\${txt}{",
"",
"\\set{cix#5}{",
"   \\@{\\w}",
"   \\format",
"      {%>6.\\` `%>5.\\` `%>5.\\` `%>10.}",
"      {{\\1}{\\3}{\\4}{\\5}}",
"   \\@{\\W\\N}",
"}",
"",
"                           }\\${html}{",
"",
"\\set{cix#5}{",
"   \\@{\\w}",
"   \\format",
"      {%*{length}{\\1}>6.\\` `%>5.\\` `%>5.\\` `%>10.}",
"      {{\\lref{\\2.html#c\\1}{\\1}}{\\3}{\\4}{\\5}}",
"   \\@{\\W\\N}",
"}",
"",
"                           }}",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   element index entry.",
"   \\\\eix{tab}{el}{cl}{fname}{nbi-nbo}{self}{cvg}{xi}{xo}",
"                            ----------------------------",
"",
"",
"                           \\formatted{",
"                           \\${txt}{",
"",
"\\set{eix#9}{",
"   \\@{\\w}",
"   \\if{\\cmp{ne}{\\1}{}}{[\\1]\\nl}{}",
"   \\format",
"      {  %>6.\\` `",
"         %>5.\\` `",
"         %@{-}{7}12.\\` `",
"         %>4.\\` `",
"         %@{-}{7}12.\\` `",
"         %@{-}{7}12.\\` `",
"         %@{-}{7}12.",
"      }",
"      { {\\2}{\\3}{\\5}{\\6}{\\7}{\\8}{\\9} }",
"   \\@{\\W}",
"}",
"",
"                           }\\${html}{",
"",
"\\set{eix#9}{",
"   \\@{\\w}",
"   \\if{\\cmp{ne}{\\1}{}}{[\\1]\\nl}{}",
"   \\format",
"      {  %*{length}{\\2}>6.\\` `",
"         %*{length}{\\3}>5.\\` `",
"         %@{-}{7}12.\\` `",
"         %>4.\\` `",
"         %@{-}{7}12.\\` `",
"         %@{-}{7}12.\\` `",
"         %@{-}{7}12.",
"      }",
"      {  {\\lref{\\4.html#e\\2}{\\2}}",
"         {\\lref{\\4.html#c\\3}{\\3}}",
"         {\\5}{\\6}{\\7}{\\8}{\\9}",
"      }",
"   \\@{\\W}",
"}",
"",
"                           }}",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"                           \\formatted{",
"",
"\\set{sec_inner#1}{",
"   \\par",
"   \\enref{in\\1}{}---\\nl\\,",
"   Inner\\` `nodes",
"   \\${html}{",
"      \\` `\\iref{ou\\1}{Outer}",
"      \\` `\\iref{c\\1}{Top}",
"      \\` `\\iref{bo\\1}{Bottom}",
"   }",
"   \\nl---\\nl",
"}",
"\\set{sec_outer#1}{",
"   \\par",
"   \\enref{ou\\1}{}---\\nl\\,",
"   Outer\\` `nodes",
"   \\${html}{",
"      \\` `\\iref{in\\1}{Inner}",
"      \\` `\\iref{c\\1}{Top}",
"      \\` `\\iref{bo\\1}{Bottom}",
"   }",
"   \\nl---\\nl",
"}",
"\\set{sec_outer_end#1}{",
"   \\${html}{",
"      \\nl",
"      \\enref{bo\\1}{}---\\nl\\,",
"      \\` `\\iref{in\\1}{Inner}",
"      \\` `\\iref{ou\\1}{Outer}",
"      \\` `\\iref{c\\1}{Top}",
"      \\` `\\iref{bo\\1}{Bottom}",
"      \\nl",
"   }",
"}",
"",
"\\set{sec_inner_end#1}{}",
"",
"                           }",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   cluster listing entry",
"   \\\\cl{cl}{sz}{coh}{cov}{prev}{prevname}{next}{nextname}",
"",
"\\${txt}{",
"\\set{cl#8}{_\\1\\,_ prev[_\\5\\,_ \\6.txt] next[_\\7\\,_ \\8.txt]\\par",
"Cluster \\1 sz \\2 self \\3 cov \\4\\nl}",
"}",
"\\${html}{\\set{cl#8}{\\",
"\\lref{\\6.html#c\\5}{<} \\lref{\\8.html#c\\7}{>}\\",
" \\enref{c\\1}{Cluster \\1} sz \\2 self \\3 cov \\4}}",
"",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   cluster cluster hits",
"   \\\\c2c{self}{{cl}{fname}{proj} ..}",
"",
"                           \\formatted{",
"                           \\${txt}{",
"",
"\\set{c2c#2}{\\@{\\w}\\format{%>6.\\` `%>4.}{{self:}{\\1}}",
"\\apply{_#3\\!{{\\format{%>6.\\` `%>4.}{{\\1:}{\\3}}\\nl}}}{\\2}\\@{\\W}}",
"",
"                           }",
"                           \\${html}{",
"\\set{c2c#2}{",
"   \\@{\\P\\w}",
"   \\format{%20@{:}{14}.}{{self:\\` `\\1}}",
"   \\nl",
"   \\apply",
"      {_#3\\!{{",
"         \\format{",
"            %20@{:}{14}*{length}{cl\\` `\\1:\\` `\\3}.",
"         }{",
"            {\\lref{\\2.html#c\\1}{cl\\` `\\1}:\\` `\\3}",
"         }",
"         \\nl",
"         }}",
"      }",
"      {\\2}",
"      \\@{\\W}",
"}",
"                           }}",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   element inner self",
"   \\\\eis{tab}{el}{nbi-nbo}{self}{cvg}{xi}{xo}{SUM}",
"                 ----------------------------",
"",
"                           \\formatted{",
"\\set{eis#8}{",
"   \\@{\\w}",
"   \\if{\\cmp{ne}{\\1}{}}{[\\1]\\nl}{}",
"   \\${html}{\\enref{e\\2}{}}\\",
"   \\format",
"      {%>6.\\` `%@{-}{7}12.\\` `%>4.\\` `%@{-}{7}12.\\` `%@{-}{7}12.\\` `%@{-}{7}12.}",
"      { {\\2}{\\3}{\\4}{\\5}{\\6}{\\7} }",
"      \\@{\\S}<\\8>",
"   \\@{\\W}",
"}",
"                           }",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"",
"   element inner other",
"   \\\\eio{el}{cl}{fname}{clsz/ct}{proj}{cvg}{xi}",
"",
"                           \\formatted{",
"                           \\${txt}{",
"\\set{eio#7}{",
"   \\@{\\w}",
"   \\format",
"      {%>8.\\` `%@{/}{5}10.\\` `%>4.\\` `%@{-}{7}12.\\` `%@{-}{7}12.}",
"      { {\\2}{\\4}{\\5}{\\6}{\\7} }",
"   \\@{\\W}",
"}",
"                           }\\${html}{",
"",
"\\set{eio#7}{",
"   \\@{\\w}",
"   \\format",
"      {  %>8*{length}{\\2}.\\` `",
"         %@{/}{5}10.\\` `%>4.\\` `",
"         %@{-}{7}12.\\` `",
"         %11.",
"      }",
"      {  {\\lref{\\3.html#c\\2}{\\2}}",
"         {\\4}{\\5}{\\6}{\\7}",
"      }",
"   \\@{\\W}",
"}",
"",
"                           }}",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   element alien self",
"   \\\\eas{tab}{el}{cl}{fname}{sz}{{nbi-nbo}{self}{cvg}{xi}{xo}}",
"                                 ----------------------------",
"",
"                           \\formatted{",
"                           \\${txt}{",
"",
"",
"\\set{eas#6}{",
"   \\@{\\w}",
"   \\if{\\cmp{ne}{\\1}{}}{[\\1]\\nl}{}",
"   \\format",
"      {%>6.\\` `%@{#}{6}12.}",
"      {{\\2}{\\3#\\5}}",
"   \\format",
"      {%@{-}{5}10.\\` `%>4.\\` `%@{-}{7}12.\\` `%@{-}{7}12.\\` `%@{-}{7}12.}",
"      {\\6}",
"   \\@{\\W}",
"}",
"",
"                           }\\${html}{",
"",
"\\set{eas#6}{",
"   \\@{\\w}",
"   \\if{\\cmp{ne}{\\1}{}}{[\\1]\\nl}{}",
"   \\format",
"      {  %*{length}{\\2}>6.\\` `",
"         %*{length}{\\3#\\5}@{#}{6}12.",
"      }",
"      {  {\\lref{\\4.html#e\\2}{\\2}}",
"         {\\lref{\\4.html#c\\3}{\\3}#\\5}",
"      }",
"   \\format",
"      {  %@{-}{5}10.\\` `",
"         %>4.\\` `",
"         %@{-}{7}12.\\` `",
"         %@{-}{7}12.\\` `",
"         %@{-}{7}12.",
"      }",
"      {\\6}",
"   \\@{\\W}",
"}",
"",
"                           }}",
"",
"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"",
"",
"   element alien this",
"   \\\\eat{ct}{proj}{cvg}{xi}{SUM}",
"",
"\\${txt}{",
"\\set{eat#5}{\\@{\\w}\\format{%22.%<7. %>4. %@{-}{7}12. %@{-}{7}12.}{",
"   {}{/\\1}{\\2}{\\3}{\\4}}\\@{\\S}<\\5>\\@{\\W}\\nl}",
"}",
"\\${html}{",
"\\set{eat#5}{\\@{\\w}\\format{%22.%<7. %>4. %@{-}{7}12. %@{-}{7}12.}{",
"   {}{/\\1}{\\2}{\\3}{\\4}}\\@{\\S}<\\5>\\@{\\W}}",
"}",
"}",
"",
"",
"",
NULL
}  ;