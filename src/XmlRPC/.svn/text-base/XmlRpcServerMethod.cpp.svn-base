#include "XmlRPC/XmlRpcServer.h"
#include "XmlRPC/XmlRpcServerMethod.h"

namespace XmlRpc
{
	XmlRpcServerMethod::XmlRpcServerMethod(std::string const& name, XmlRpcServer* server)
	{
		_name = name;
		_server = server;
		if (_server)
		{
			_server->addMethod(this);
		}
	}
	XmlRpcServerMethod::~XmlRpcServerMethod()
	{
		if (_server)
		{
			_server->removeMethod(this);
		}
	}
} // namespace XmlRpc
