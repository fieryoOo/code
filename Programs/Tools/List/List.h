#include <cstdio>
#include <cstdlib>

template <class T>
class ListBase {
public:
   virtual ~ListBase() {}

//   virtual void flush() = 0;

//   virtual void Push(const T data) = 0;
   virtual void Append(const T data) = 0;
//   virtual void DelFront() = 0;

//   virtual const T GetFirst() const = 0;
   virtual T GetFirst() = 0;
//   virtual const T GetLast() const = 0;
   virtual T GetLast() = 0;

//   virtual const int IsEmpty() const = 0;
};

template <class T>
class DataNode {
   T data;
public:
   DataNode *Last, *Next;
   DataNode(T datain) {
      data = datain;
      Last = Next = NULL;
   }
   T Data() {
      return data;
   }
};

template <class T>
class List : ListBase<T> {
   DataNode<T> *head, *tail;

public:
   List() {}

   void Append(const T data) {
      tail = new DataNode<T>(data);
      if( ! head ) head = tail;
   }

   T GetFirst() {
      if( ! head ) return NULL;
      return (*head).Data();
   }

   T GetLast() {
      if( ! tail ) return NULL;
      return (*tail).Data();
   }

   ~List() {}
};


/*
template <class T> void
List<T>::Append(const T data) {
   tail = new DataNode<T>(data);
   if( ! head ) head = tail;
};

template <class T> T
List<T>::GetFirst() {
   if( ! head ) return NULL;
   return (*head).Data();
};
*/
