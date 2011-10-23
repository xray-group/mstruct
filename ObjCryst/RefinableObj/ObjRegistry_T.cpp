
#include "RefinableObj.h"

namespace ObjCryst {

//######################################################################
//    ObjRegistry
//######################################################################
#ifdef __WX__CRYST__
bool operator==(const wxString&wx,const string&str)
{
   return wx==str.c_str();
}
bool operator==(const string&str,const wxString&wx)
{
   return wx==str.c_str();
}
#endif

template<class T> ObjRegistry<T>::ObjRegistry():
mName("")
#ifdef __WX__CRYST__
,mpWXRegistry(0)
#endif
{
   VFN_DEBUG_MESSAGE("ObjRegistry::ObjRegistry()",5)
}

template<class T> ObjRegistry<T>::ObjRegistry(const string &name):
mName(name)
#ifdef __WX__CRYST__
,mpWXRegistry(0)
#endif
{
   VFN_DEBUG_MESSAGE("ObjRegistry::ObjRegistry(name):"<<mName,5)
}

//:TODO: a copy constructor
template<class T> ObjRegistry<T>::~ObjRegistry()
{
   VFN_DEBUG_MESSAGE("ObjRegistry::~ObjRegistry():"<<mName,5)
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
}

template<class T> void ObjRegistry<T>::Register(T &obj)
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::Register():"<<obj.GetName(),2)
   typename vector<T*>::iterator pos=find(mvpRegistry.begin(),mvpRegistry.end(),&obj);
   if(pos!=mvpRegistry.end())
   {
      VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Register():"<<obj.GetName()<<"Already registered!",2)
      return;
   }
   mvpRegistry.push_back(&obj);
   mListClock.Click();
   #ifdef __WX__CRYST__
   if(0!=mpWXRegistry) 
      mpWXRegistry->Add(obj.WXCreate(mpWXRegistry));
   #endif
   //this->Print();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Register():"<<obj.GetName(),2)
}

template<class T> void ObjRegistry<T>::DeRegister(T &obj)
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::Deregister(&obj)",2)
   //this->Print();
   typename vector<T*>::iterator pos=find(mvpRegistry.begin(),mvpRegistry.end(),&obj);
   if(pos==mvpRegistry.end())
   {
      VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(&obj):NOT FOUND !!!",2)
      return; //:TODO: throw something ?
   }
   #ifdef __WX__CRYST__
   if(0!=mpWXRegistry) mpWXRegistry->Remove(obj.WXGet());
   #endif
   mvpRegistry.erase(pos);
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(&obj)",2)
}

template<class T> void ObjRegistry<T>::DeRegister(const string &objName)
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::Deregister(name):"<<objName,2)
   
   const long i=this->Find(objName);
   if(-1==i)
   {
      VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(name): NOT FOUND !!!",2)
      return; //:TODO: throw something ?
   }
   //:KLUDGE: should directly do an iterator search on the name...
   typename vector<T*>::iterator pos=find(mvpRegistry.begin(),mvpRegistry.end(),mvpRegistry[i]);
   
   #ifdef __WX__CRYST__
   if(0!=mpWXRegistry) mpWXRegistry->Remove((*pos)->WXGet());
   #endif
   mvpRegistry.erase(pos);
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(name):",2)
}

template<class T> void ObjRegistry<T>::DeRegisterAll()
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::DeRegisterAll():",5)
   #ifdef __WX__CRYST__
   if(0!=mpWXRegistry)
   {
      typename vector<T*>::iterator pos;
      for(pos=mvpRegistry.begin();pos!=mvpRegistry.end();++pos)
         mpWXRegistry->Remove((*pos)->WXGet());
   }
   #endif
   mvpRegistry.clear();
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::DeRegisterAll():",5)
}

template<class T> void ObjRegistry<T>::DeleteAll()
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::DeleteAll():",5)
   vector<T*> reg=mvpRegistry;//mvpRegistry will be modified as objects are deleted, so use a copy
   typename vector<T*>::iterator pos;
   for(pos=reg.begin();pos!=reg.end();++pos) delete *pos;
   mvpRegistry.clear();
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::DeleteAll():",5)
}

template<class T> T& ObjRegistry<T>::GetObj(const unsigned int i)
{
   return *(mvpRegistry[i]);
}

template<class T> const T& ObjRegistry<T>::GetObj(const unsigned int i) const
{
   return *(mvpRegistry[i]);
}

template<class T> T& ObjRegistry<T>::GetObj(const string &objName)
{
   const long i=this->Find(objName);
   return *(mvpRegistry[i]);
}

template<class T> const T& ObjRegistry<T>::GetObj(const string &objName) const
{
   const long i=this->Find(objName);
   return *(mvpRegistry[i]);
}

template<class T> T& ObjRegistry<T>::GetObj(const string &objName,
                                                  const string& className)
{
   const long i=this->Find(objName,className);
   return *(mvpRegistry[i]);
}

template<class T> const T& ObjRegistry<T>::GetObj(const string &objName,
                                                        const string& className) const
{
   const long i=this->Find(objName,className);
   return *(mvpRegistry[i]);
}

template<class T> long ObjRegistry<T>::GetNb()const{return (long)mvpRegistry.size();}

template<class T> void ObjRegistry<T>::Print()const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Print():",2)
   cout <<mName<<" :"<<this->GetNb()<<" object registered:" <<endl;
   
   for(long i=0;i<this->GetNb();++i)
      cout <<i<<"("<<this->GetObj(i).GetName()<<")"<<endl;
}

template<class T> void ObjRegistry<T>::SetName(const string &name){ mName=name;}

template<class T> const string& ObjRegistry<T>::GetName()const { return mName;}

template<class T> long ObjRegistry<T>::Find(const string &objName) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(objName)",2)
   long index=-1;
   //bool error=false;
   for(long i=this->GetNb()-1;i>=0;i--) 
      if( mvpRegistry[i]->GetName() == objName) return i;
   //      if(-1 != index) error=true ;else index=i;
   //if(true == error)
   //{
   //   cout << "ObjRegistry::Find(name) : ";
   //   cout << "found duplicate name ! This *cannot* be !!" ;
   //   cout << objName <<endl;
   //   this->Print();
   //   throw 0;
   //}
   cout << "ObjRegistry<T>::Find("<<objName<<"): Not found !!"<<endl;
   this->Print();
   throw ObjCrystException("ObjRegistry<T>::Find("+objName+"): Not found !!");
   return index;
}

template<class T> long ObjRegistry<T>::Find(const string &objName,
                                            const string &className,
                                             const bool nothrow) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(objName,className)",2)
   long index=-1;
   //bool error=false;
   for(long i=this->GetNb()-1;i>=0;i--) 
      if( mvpRegistry[i]->GetName() == objName) 
         if(className==mvpRegistry[i]->GetClassName()) return i;
   //      if(-1 != index) error=true ;else index=i;
   //if(true == error)
   //{
   //   cout << "ObjRegistry::Find(name) : ";
   //   cout << "found duplicate name ! This *cannot* be !!" ;
   //   cout << objName <<endl;
   //   this->Print();
   //   throw 0;
   //}
   cout << "ObjRegistry<T>::Find("<<objName<<","<<className<<"): Not found !!"<<endl;
   this->Print();
   if(nothrow==false)
      throw ObjCrystException("ObjRegistry<T>::Find("+objName+","+className+"): Not found !!");
   return index;
}

template<class T> long ObjRegistry<T>::Find(const T &obj) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(&obj)",2)
   for(long i=this->GetNb()-1;i>=0;i--) 
      if( mvpRegistry[i]== &obj)  return i;
   //:TODO: throw something
   return -1;
}

template<class T> long ObjRegistry<T>::Find(const T *pobj) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(&obj)",2)
   for(long i=this->GetNb()-1;i>=0;i--) 
      if( mvpRegistry[i]== pobj)  return i;
   //:TODO: throw something
   return -1;
}

template<class T> const RefinableObjClock& ObjRegistry<T>::GetRegistryClock()const{return mListClock;}

#ifdef __WX__CRYST__
template<class T> WXRegistry<T>* ObjRegistry<T>::WXCreate(wxWindow *parent)
{
   VFN_DEBUG_MESSAGE("ObjRegistry<T>::WXCreate()",2)
   mpWXRegistry=new WXRegistry<T> (parent,this);
   for(int i=0;i<this->GetNb();i++) 
      mpWXRegistry->Add(this->GetObj(i).WXCreate(mpWXRegistry));
   return mpWXRegistry;
}
template<class T> void ObjRegistry<T>::WXDelete()
{
   if(0!=mpWXRegistry)
   {
      VFN_DEBUG_MESSAGE("ObjRegistry<T>::WXDelete()",2)
      delete mpWXRegistry;
   }
   mpWXRegistry=0;
}
template<class T> void ObjRegistry<T>::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("ObjRegistry<T>::WXNotifyDelete()",2)
   mpWXRegistry=0;
}
#endif

}

