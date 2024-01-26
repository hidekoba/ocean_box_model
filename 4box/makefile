#生成するプログラム名（複数可）
PROGS=a.out

#使用するソース
SRC=chemeq.f90 o2.f90 co2.f90 main.f90

#ソース名前→オブジェクト
OBJS=$(SRC:.f90=.o)

#コンパイラ名
FC=gfortran

#コンパイルオプション
FCFLAGS=-Wall -O2

#サフィックス
.SUFFIXES: .f90

all:$(PROGS)

a.out:$(OBJS)
	$(FC) $(FCFLAGS) -o $@ $(OBJS)

#main の依存
co2.o:chemeq.o
main.o:chemeq.o co2.o

.f90.o:
	$(FC) $(FCFLAGS) -c $< $(FCFLAGS) -o $@

clean:
	rm -f $(OBJS) $(OBJS:.o=.mod) $(PROGS)
