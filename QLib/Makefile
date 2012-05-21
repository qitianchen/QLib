INC = ..\eigen3
CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -I$(INC)

OBJS =		main.o TwoFactorsBlackScholesMertonEquation.o ExchangePutPayoff.o \
ExchangePutBoundaries.o BoundaryCall.o BarrierBoundaries.o \
BlackScholesMertonEquation.o BarrierPayoff.o \
PayOffCall.o

LIBS =

TARGET =	QLib_Demo.exe

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
