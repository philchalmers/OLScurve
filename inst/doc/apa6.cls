% LaTeX class for APA 6 manuscripts and document styled drafts
% Author: Phil Chalmers
% Lasted edited on Nov 20, 2011

%Options to the class:
    % man - manuscript mode (default)
    % doc - document mode
    % noendfloat - inline figures and tables
    % timesnr - times new roman font

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{apa6}[2011/11/21 LaTeX class for APA 6 manuscripts]
\RequirePackage{ifthen}

%Constants
\def\def@noendfloat{false}
\def\def@timesnr{false}

%Options to class

\DeclareOption{noendfloat}{
    \def\def@noendfloat{true}
}

\DeclareOption{timesnr}{
    \def\def@timesnr{true}
}

\DeclareOption{man}{
    \def\def@doc{manmode}
}

\DeclareOption{doc}{
    \def\def@doc{docmode}
}

\ProcessOptions


\@ifundefined{def@doc}{
     \def\def@doc{manmode}
     \ClassWarningNoLine{apa6}{Using default mode (man)}
}{}

% Change greek symbols to upper case

\RequirePackage{upgreek} %makes greek letter not in italics, use \up____
\renewcommand{\alpha}     {\upalpha     }
\renewcommand{\beta}      {\upbeta      }
\renewcommand{\gamma}     {\upgamma     }
\renewcommand{\delta}     {\updelta     }
\renewcommand{\epsilon}   {\upepsilon   }
\renewcommand{\varepsilon}{\upvarepsilon}
\renewcommand{\zeta}      {\upzeta      }
\renewcommand{\eta}       {\upeta       }
\renewcommand{\theta}     {\uptheta     }
\renewcommand{\vartheta}  {\upvartheta  }
\renewcommand{\iota}      {\upiota      }
\renewcommand{\kappa}     {\upkappa     }
\renewcommand{\lambda}    {\uplambda    }
\renewcommand{\mu}        {\upmu        }
\renewcommand{\nu}        {\upnu        }
\renewcommand{\xi}        {\upxi        }
\renewcommand{\pi}        {\uppi        }
\renewcommand{\varpi}     {\upvarpi     }
\renewcommand{\rho}       {\uprho       }
\renewcommand{\varrho}    {\upvarrho    }
\renewcommand{\sigma}     {\upsigma     }
\renewcommand{\varsigma}  {\upvarsigma  }
\renewcommand{\tau}       {\uptau       }
\renewcommand{\upsilon}   {\upupsilon   }
\renewcommand{\phi}       {\upphi       }
\renewcommand{\varphi}    {\upvarphi    }
\renewcommand{\chi}       {\upchi       }
\renewcommand{\psi}       {\uppsi       }
\renewcommand{\omega}     {\upomega     }

%man
\ifthenelse{\equal{\def@doc}{manmode}}{
    \LoadClass[12pt, \CurrentOption]{article}
    \RequirePackage[T1]{fontenc}
    \RequirePackage[american]{babel}
    \RequirePackage{csquotes}
    \RequirePackage[style=apa]{biblatex}
    \DeclareLanguageMapping{american}{american-apa}
    \RequirePackage{url}
    \RequirePackage{graphicx}
    \RequirePackage{amssymb}
    \RequirePackage{amsmath}
    \RequirePackage{booktabs} %use \midrule only (not \toprule \bottomrule)
    \RequirePackage{rotating}
    \RequirePackage{verse}
    \RequirePackage{enumitem}
    \setlist{nolistsep}
    \ifthenelse{\equal{\def@noendfloat}{true}}{}{
        \RequirePackage[nofiglist,notablist,tablesfirst, nofighead, notabhead]{endfloat}
        \renewcommand{\figureplace}{\begin{center} ---INSERT FIGURE~\thepostfig\
            HERE--- \end{center}}
        \renewcommand{\tableplace}{\begin{center} ---INSERT TABLE~\theposttbl\
            HERE--- \end{center}}
    }
    \RequirePackage{endnotes}
    \let\footnote=\endnote
    \def\enotesize{\normalsize}
    \RequirePackage{indentfirst}
    \RequirePackage[margin=1in,top=1.25in]{geometry}
    \RequirePackage{fancyhdr}
    \RequirePackage{setspace}
    \RequirePackage[config=mycaption]{caption}
    \DeclareCaptionLabelSeparator{period}{. }
    \DeclareCaptionLabelSeparator{newline}{\newline}
    \captionsetup[table]{singlelinecheck=false,labelfont=default,labelsep=newline, textfont=it,justification={raggedright,doublespacing}}
    \captionsetup[figure]{style=default,labelsep=period, labelfont=it,justification={raggedright,doublespacing}}
    \setlength{\parskip}{0ex plus 0ex minus 0ex}
    \doublespacing
    \pagestyle{fancy}
    \renewcommand{\headrulewidth}{0pt}
    \setlength{\parindent}{.5in}
    \fancyhead[R]{}
    \fancyfoot{}
    \makeatletter
    \renewcommand\section{\@startsection{section}{1}{0mm}{\lineskip}{.2\baselineskip}{\normalfont\bfseries\centering}}
    \renewcommand\subsection{\@startsection{subsection}{1}{0mm}{\lineskip}{.1\baselineskip}{\normalfont\bfseries\raggedright}}
    \renewcommand\subsubsection[1]{\@startsection{subsubsection}{1}{.5in}{\lineskip}{.0\baselineskip}{\normalfont\bfseries\raggedright}{#1. }}
    \makeatother
    \setcounter{secnumdepth}{0}
    \textwidth = 450pt
    \oddsidemargin = 10pt
    \newcommand{\vlength}[0]{449pt}
    \raggedright
    \def\runningheadvar{RUNNING HEADER}
    \newcommand{\runninghead}[1]{\def\runningheadvar{#1}}
    \renewcommand{\title}[1]{\def\titlevar{#1}}
    \renewcommand{\author}[1]{\def\authorvar{#1}}
    \def\authornotevar{}
    \newcommand{\authornote}[1]{\def\authornotevar{\section{Author Note} \vspace{.3cm} \indent #1}}
    \def\abstractvar{}
    \def\keywordsvar{}
    \def\authorinfovar{}
    \renewcommand{\abstract}[1]{\def\abstractvar{\section{Abstract}\noindent #1 \newline}}
    \newcommand{\keywords}[1]{\def\keywordsvar{\indent \emph{Keywords}: #1}}
    \newcommand{\authorinfo}[1]{\def\authorinfovar{#1}}
    \setlength{\parindent}{.5in}
    \renewcommand{\maketitle}{
        \fancyhead[L]{Running Head: \MakeUppercase\runningheadvar\hfill\thepage}
        \vspace*{3cm}
        \begin{center}
        \titlevar \\ \authorvar \\ \authorinfovar
        \end{center}
        \vspace{4cm}
        \authornotevar
        \newpage
        \fancyhead[L]{\MakeUppercase\runningheadvar\hfill\thepage}
        \abstractvar \keywordsvar
        \newpage
        \begin{center}
        \titlevar
        \end{center}
    }
    \ifthenelse{\equal{}{}}{\RequirePackage{txfonts}}{}
}{}

%doc
\ifthenelse{\equal{\def@doc}{docmode}}{
    \LoadClass[\CurrentOption]{article}
    \usepackage[margin=1.5in]{geometry}
    \usepackage[T1]{fontenc}
    \usepackage[american]{babel}
    \usepackage{csquotes}
    \usepackage[style=apa]{biblatex}
    \DeclareLanguageMapping{american}{american-apa}    
    \usepackage{url}
    \usepackage{graphicx}
    \usepackage{amssymb,amsmath}
    \usepackage{booktabs} %use \midrule only (not \toprule \bottomrule)
    \usepackage{rotating}
    \usepackage{verse}
    \usepackage{setspace}
    \usepackage[config=mycaption]{caption}
    \usepackage{fancyhdr}
    \usepackage{ifthen}
    \usepackage{enumitem}
    \setlist{nolistsep}
    \pagestyle{fancy}
    \DeclareCaptionLabelSeparator{period}{. }
    \DeclareCaptionLabelSeparator{newline}{\newline}
    \captionsetup[table]{singlelinecheck=false,labelfont=default,labelsep=newline, textfont=it,justification={raggedright,doublespacing}}
    \captionsetup[figure]{style=default,labelsep=period, labelfont=it,justification={raggedright,doublespacing}}
    \setlength{\parskip}{0ex plus 0ex minus 0ex}
    \setlength{\parindent}{.35in}
    \makeatletter
    \newcommand{\theendnotes}[0]{}
    \def\authornotevar{}
    \newcommand{\authornote}[1]{\def\authornotevar{#1}}
    \renewcommand\section{\@startsection {section}{1}{\z@}%
      {-2.5ex \@plus -1ex \@minus -.2ex}%
      {2ex \@plus .2ex}%
      {\centering\normalfont\large\bfseries}}
    \makeatother
    \setcounter{secnumdepth}{0}
    \newcommand{\vlength}[0]{349pt}
    \def\runningheadvar{RUNNING HEAD}
    \def\titlevar{Title}
    \def\authorvar{\normalsize Authors \\ Departments \\ Institutions}
    \def\abstractvar{}
    \def\keywordsvar{}
    \newcommand{\runninghead}[1]{\def\runningheadvar{#1}}
    \renewcommand{\title}[1]{\def\titlevar{\Large #1}}
    \renewcommand{\author}[1]{\def\authorvar{\normalsize #1}}
    \renewcommand{\abstract}[1]{\def\abstractvar{\section{Abstract} \noindent #1}}
    \newcommand{\keywords}[1]{\def\keywordsvar{\newline \indent \emph{Keywords}: #1}}
    \newcommand{\authorinfo}[1]{\def\authorinfovar{#1}}
    \fancyhead{}
    \fancyfoot{}
    \renewcommand{\headrulewidth}{0pt}
    \renewcommand{\footrulewidth}{0pt}
    \fancyhead[C]{\MakeUppercase\runningheadvar}
    \fancyfoot[C]{\thepage}
    \renewcommand{\maketitle}{
        \thispagestyle{plain}
        \begin{center} \titlevar \\ \medskip \authorvar \\ \medskip Last compiled on \today \end{center}
        \let\thefootnote\relax\footnotetext{\authornotevar}
        \abstractvar \keywordsvar
        \newpage
        \section{Introduction}
    }
}{}

%%%Definitions after packages loaded for all styles%%%
% \pkg{} and \proglang{} commands
\let\proglang=\textsf
\newcommand{\pkg}[1]{{\fontseries{b}\selectfont #1}}

% citation style
\let\citeA=\textcite
\let\cite=\parencite
\let\citep=\parencite
\let\citeNP=\nptextcite
\newcommand{\citeAs}[1]{\citeauthor{#1}'s (\citeyear{#1})}
\newcommand{\citeyearpar}[1]{(\citeyear{#1})}
